#! /usr/bin/env python3
"""Keeps one transcript per gene of a GFF3 file, the one with the longest spliced CDS (the first in the
file on ties), with its gene line and its exon, CDS and UTR lines, e.g. as the reference for comparing a
gene prediction with gffcompare. With --no-utrs, the exon and UTR lines are left out as well, so every
transcript keeps only its CDS lines.

The lines are grouped into genes and transcripts as GeenuFF's import groups them; lines it leaves out
(see the import summary) are left out here too, as are non-coding genes and every line of other types.
The kept lines are written as they are, in file order.

With --add-stop-codons, a stop codon missing from a CDS is added as the import adds it, before the
longest transcript is picked (same rules as the importer). The last CDS line is extended, or a CDS
line added past an intron after it; where the codon lies past the exons, the last exon line, the
transcript line and, if need be, the gene line are extended too. A UTR line holding the codon is
shortened to start after it, or left out where the codon covers it fully, and counted.
"""
import argparse
import logging
import os
from collections import Counter
from collections.abc import Iterator
from dataclasses import dataclass, field

from geenuff.applications.import_summary import DROPPED_LINE_REASONS
from geenuff.applications.importer import (GENE_TYPES, PIECE_KEYS, TRANSCRIPT_TYPES, FastaImporter,
                                           MissingGeneLines, OrganizedGeenuffImporterGroup, OrganizedGFFEntries)
from geenuff.base import orm
from geenuff.base.helpers import (Compression, NonMatchableIDs, add_output_compression_args, check_output_paths,
                                  insert_before_extension, match_seqids, merge_gff_spans, open_out, output_compression,
                                  output_path, require_fasta_and_gff, require_input_files, setup_logging,
                                  spliced_cds_sequence, strand_or_none)
from geenuff.parsing.fastahelper import FastaParser
from geenuff.parsing.gffhelper import GFFObject, GFFParser

logger = logging.getLogger(__name__)

# inserted before the extension of the GFF3 file's name for the default output
OUT_SUFFIX = '_longest'
NO_UTRS_SUFFIX = '_longest_cds_only'

# {(kind: 'gene', 'transcript', 'exons', 'cds' or 'utrs', ID of the gene or transcript or of the transcript
# the line names, old start, old end): (new start, new end), None for a line left out}
Respans = dict[tuple[str, str, int, int], tuple[int, int] | None]
# {(transcript ID, start, end of the CDS line they follow): [(start, end, phase) of the CDS lines added]}
AddedCds = dict[tuple[str, int, int], list[tuple[int, int, int]]]


@dataclass
class Selection:
    """The genes and transcripts kept and the changes to their lines.

    Attributes:
        genes: IDs of the genes kept
        transcripts: {ID: sequence ID} of their longest coding transcripts
        respans: new spans of the lines changed by adding stop codons
        added: CDS lines added past an intron by adding stop codons
        outcomes: genes per outcome, written, non_coding or no_strand (see the module docstring), and
            written transcripts with a stop codon added (stop_codon, stop_codon_past_exons)
        utr_lines_left_out: UTR lines of written transcripts left out for being covered by the stop codon
        stop_codon_in_utr: IDs of the written transcripts with a UTR line holding the stop codon added
    """
    genes: set[str] = field(default_factory=set)
    transcripts: dict[str, str] = field(default_factory=dict)
    respans: Respans = field(default_factory=dict)
    added: AddedCds = field(default_factory=dict)
    outcomes: Counter = field(default_factory=Counter)
    utr_lines_left_out: int = 0
    stop_codon_in_utr: list[str] = field(default_factory=list)


def coding_length(cds: list[GFFObject], exons: list[GFFObject]) -> int:
    """Bases from the first to the last CDS line that lie in the exons: the spliced CDS length the import
    picks the longest transcript by (OrganizedGeenuffImporterGroup._coding_length), here from the GFF
    lines."""
    start, end = min(c.start for c in cds), max(c.end for c in cds)
    return sum(max(0, min(end, e_end) - max(start, e_start) + 1)
               for e_start, e_end in merge_gff_spans((e.start, e.end) for e in exons))


def add_stop_codon(t: GFFObject, entries: dict, is_plus_strand: bool, coord: orm.Coordinate) -> int | None:
    """Adds a missing stop codon to a transcript's lines in place, as the import adds it (see
    OrganizedGeenuffImporterGroup.add_missing_stop_codon).

    Returns:
        how many of the 3 bases lie past the exons, None where no stop codon was added
    """
    # a transcript reaching past its sequence is masked whole by the import, without further checks
    # so it is ignored here as well (also if it's past the 3' end, nothing could be added)
    if any(line.start < 1 or line.end > coord.length for line in [t] + entries['exons'] + entries['cds']):
        return None
    cds_seq = spliced_cds_sequence(coord, entries['cds'], is_plus_strand)
    added = OrganizedGeenuffImporterGroup.add_missing_stop_codon(t, entries, is_plus_strand, coord, cds_seq)
    return None if added is None else added[1]


def record_stop_codon(t_id: str, old_cds: list[GFFObject], old_exons: list[GFFObject], entries: dict,
                      is_plus_strand: bool, respans: Respans, added: AddedCds) -> None:
    """Records the CDS and exon lines a stop codon was added to with their new spans, the CDS lines added
    past an intron, after the 3' CDS line, with their phase, and the UTR lines holding the codon, shortened
    to start after it or left out where it covers them."""
    n = len(old_cds)
    cds = entries['cds']
    # the lines added lie 3' of the old ones: after them on the plus strand, before them on the minus one
    kept, new = (cds[:n], cds[n:]) if is_plus_strand else (cds[len(cds) - n:], cds[:len(cds) - n])
    # the positions of the added codon: those of the added lines, and below those by which old lines grew
    codon = [p for line in new for p in range(line.start, line.end + 1)]
    for old, line in zip(old_cds, kept):
        # span needs to be changed on a line that was in the file beforehand
        if (old.start, old.end) != (line.start, line.end):
            respans['cds', t_id, old.start, old.end] = (line.start, line.end)
            codon += range(old.end + 1, line.end + 1) if is_plus_strand else range(line.start, old.start)
    for utr in entries['utrs']:
        # adapt UTR spans if the stop codon was at its start
        in_utr = [p for p in codon if utr.start <= p <= utr.end]
        if in_utr:
            span = (max(in_utr) + 1, utr.end) if is_plus_strand else (utr.start, min(in_utr) - 1)
            respans['utrs', t_id, utr.start, utr.end] = span if span[0] <= span[1] else None
    for old, line in zip(old_exons, entries['exons']):
        # adapt the exon spans
        if (old.start, old.end) != (line.start, line.end):
            respans['exons', t_id, old.start, old.end] = (line.start, line.end)
    if new:
        # new CDS line is necessary
        last = old_cds[-1 if is_plus_strand else 0]
        in_reading_order = new if is_plus_strand else new[::-1]
        # the bases of the stop codon from a line on, i.e. those of an earlier line skipped before a
        # codon starts in it
        lengths = [line.end - line.start + 1 for line in in_reading_order]
        added[t_id, last.start, last.end] = [(line.start, line.end, sum(lengths[i:]) % 3)
                                             for i, line in enumerate(in_reading_order)]


def sequences(fasta_path: str, gff_path: str, seqids: set[str]) -> Iterator[tuple[str, orm.Coordinate]]:
    """Yields (GFF sequence ID, Coordinate holding its sequence) for the given GFF sequence IDs, one
    FASTA sequence at a time, paired with the FASTA IDs as the import pairs them (see
    helpers.match_seqids).

    Raises:
        SystemExit: where a GFF sequence ID matches no FASTA sequence or several
    """
    try:
        pairs = match_seqids(FastaParser.ids(fasta_path), GFFParser.seqids(gff_path))
    except NonMatchableIDs as e:
        raise SystemExit(f'{e}; no file written')
    gff_ids = {fasta_id: gff_id for gff_id, fasta_id in pairs.items() if gff_id in seqids}
    for fasta_id, sequence in FastaImporter.parse_fasta(fasta_path):
        if fasta_id in gff_ids:
            yield gff_ids[fasta_id], orm.Coordinate(sequence=sequence, length=len(sequence), seqid=fasta_id)


def longest_transcripts(gff_path: str, fasta_path: str | None) -> Selection:
    """The longest coding transcript of each gene, with stop codons added where fasta_path is given.

    Raises:
        SystemExit: where coding transcripts name no Parent gene (see importer.MissingGeneLines)
    """
    organizer = OrganizedGFFEntries(gff_path)
    try:
        organizer.load_organized_entries()
    except MissingGeneLines as e:
        raise SystemExit(f'{e} No file written.')
    # the lines the import's grouping leaves out, with every line below them, per reason
    for reason, text in DROPPED_LINE_REASONS.items():
        n = (sum(organizer.stats.dropped_lines.get(reason, {}).values())
             + sum(organizer.stats.dropped_lines_below.get(reason, {}).values()))
        if n:
            logger.info(f'Left out {n} lines {text}')

    entries = organizer.organized_entries
    # with a FASTA file, the genes are handled sequence by sequence alongside it, one sequence in memory
    # at a time; without one, no sequence is needed
    coords = (sequences(fasta_path, gff_path, set(entries)) if fasta_path
              else ((seqid, None) for seqid in entries))
    selection = Selection()
    outcomes = selection.outcomes
    for seqid, coord in coords:
        for group in entries[seqid]:
            gene = group['super_locus']
            # transcripts made for CDS lines without a transcript line are no transcripts of the file
            # (side effect of using the importer's functions, floating CDS get treated by it)
            coding = [(t, t_entries) for t, t_entries in group['transcripts'].items()
                      if t_entries['cds'] and not t_entries['floating']]
            if not coding:
                outcomes['non_coding'] += 1
                continue
            # stop codons are added to every coding transcript before the longest is picked, as the
            # import does; the changes are recorded per transcript and kept only for the one written
            t_respans, t_added, n_past = {}, {}, {}
            if coord is not None and gene.strand is not None:
                for t, t_entries in coding:
                    # a transcript on neither strand is taken to be on its gene's, as the import does
                    is_plus_strand = strand_or_none(t if t.strand is not None else gene)
                    old_cds, old_exons = list(t_entries['cds']), list(t_entries['exons'])
                    old_t_span = (t.start, t.end)
                    n_past[t] = add_stop_codon(t, t_entries, is_plus_strand, coord)
                    if n_past[t] is not None:
                        record_stop_codon(t.get_ID(), old_cds, old_exons, t_entries, is_plus_strand,
                                          t_respans.setdefault(t, {}), t_added.setdefault(t, {}))
                        if (t.start, t.end) != old_t_span:
                            t_respans[t][('transcript', t.get_ID(), *old_t_span)] = (t.start, t.end)
            # max keeps the first of equal ones, as the import does
            t, t_entries = max(coding, key=lambda pair: coding_length(pair[1]['cds'], pair[1]['exons']))
            # a gene not on one definite strand, e.g. trans-spliced, is excluded from exports by the import
            lines = [gene, t] + t_entries['exons'] + t_entries['cds'] + t_entries['utrs']
            if gene.strand is None or any(line.strand != gene.strand for line in lines):
                outcomes['no_strand'] += 1
                continue
            selection.genes.add(gene.get_ID())
            selection.transcripts[t.get_ID()] = t.seqid
            outcomes['written'] += 1
            if n_past.get(t) is not None:
                selection.respans |= t_respans[t]
                selection.added |= t_added[t]
                outcomes['stop_codon'] += 1
                outcomes['stop_codon_past_exons'] += n_past[t] > 0
                utr_spans = [span for key, span in t_respans[t].items() if key[0] == 'utrs']
                if utr_spans:
                    selection.stop_codon_in_utr.append(t.get_ID())
                    selection.utr_lines_left_out += utr_spans.count(None)
                # the gene line grows with a transcript extended past it for a codon past the exons
                if t.start < gene.start or t.end > gene.end:
                    selection.respans['gene', gene.get_ID(), gene.start, gene.end] = (min(gene.start, t.start),
                                                                                      max(gene.end, t.end))
    return selection


def write_longest(gff_path: str, out_path: str, selection: Selection, no_utrs: bool,
                  compression: Compression | None) -> None:
    """Writes the comment lines and the lines of the kept genes and transcripts, exon and UTR lines only
    without no_utrs, with the spans, added CDS lines and UTR lines left out of the stop codons added; an
    exon, CDS or UTR line names only the kept transcript as Parent."""
    genes, transcripts, respans = selection.genes, selection.transcripts, selection.respans
    with open_out(out_path, compression) as f_out:
        for line in GFFParser.lines(gff_path):
            columns = line.rstrip('\n').split('\t')
            if len(columns) < 9 or line.startswith('#'):
                f_out.write(line)
                continue
            feature_type, attributes = columns[2], GFFObject.parse_attributes(columns[8])
            if feature_type in GENE_TYPES:
                kind, owner = 'gene', attributes.get('ID')
                keep = owner in genes
            elif feature_type in TRANSCRIPT_TYPES:
                kind, owner = 'transcript', attributes.get('ID')
                keep = owner in transcripts
            elif feature_type in PIECE_KEYS and not (no_utrs and PIECE_KEYS[feature_type] != 'cds'):
                parents = GFFObject.parents(attributes)
                # a line on another sequence than its transcript is left out, as the import does
                kept = [parent for parent in parents if transcripts.get(parent) == columns[0]]
                keep = bool(kept)
                # a line shared with isoforms not kept names only the kept transcript
                if keep and len(kept) < len(parents):
                    columns[8] = GFFObject.with_parents(columns[8], kept)
                    line = '\t'.join(columns) + '\n'
                kind, owner = PIECE_KEYS[feature_type], kept[0] if kept else None
            else:
                keep = False
            if not keep:
                continue
            # a line changed by an added stop codon is found by its kind, its gene or transcript and its span
            # in the file
            key = (kind, owner, int(columns[3]), int(columns[4]))
            if key in respans:
                if respans[key] is None:
                    continue
                columns[3:5] = [str(position) for position in respans[key]]
                line = '\t'.join(columns) + '\n'
            f_out.write(line)
            # CDS lines added past an intron follow the 3' CDS line, copying all its columns but span and phase
            for start, end, phase in selection.added.get(key[1:], []) if kind == 'cds' else []:
                f_out.write('\t'.join(columns[:3] + [str(start), str(end)] + columns[5:7] + [str(phase), columns[8]])
                            + '\n')


def main(args: argparse.Namespace) -> None:
    setup_logging(args.log_file)
    logger.info(f'Keeping the longest coding transcript per gene of "{os.path.basename(args.gff)}"'
                f'{", without exon and UTR lines" if args.no_utrs else ""}'
                f'{", stop codons added from " + os.path.basename(args.fasta) if args.add_stop_codons else ""}')
    selection = longest_transcripts(args.gff, args.fasta if args.add_stop_codons else None)
    outcomes = selection.outcomes
    if outcomes['non_coding']:
        logger.info(f'Left out {outcomes["non_coding"]} genes without a coding transcript')
    if outcomes['no_strand']:
        logger.info(f'Left out {outcomes["no_strand"]} coding genes whose lines are not all on one definite strand')
    if args.add_stop_codons:
        logger.info(f'Added a stop codon to {outcomes["stop_codon"]} of the transcripts written, '
                    f'{outcomes["stop_codon_past_exons"]} of them past their exons')
    in_utr = selection.stop_codon_in_utr
    if in_utr:
        logger.warning(f'{len(in_utr)} of them have a UTR line holding the stop codon, shortened to start after it '
                       f'({selection.utr_lines_left_out} UTR lines covered by it left out), e.g. '
                       f'{", ".join(in_utr[:5])}; a UTR annotated to start with a stop codon likely comes from a '
                       f'conversion that put the stop codon into the UTR instead of the CDS')
    write_longest(args.gff, args.out, selection, args.no_utrs, output_compression(args, args.gff))
    logger.info(f'Wrote "{args.out}": {outcomes["written"]} genes, one transcript each')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--gff', type=str, required=True, help='GFF3 file (may be gzipped)')
    parser.add_argument('--fasta', type=str, help='FASTA file of the genome (may be gzipped), for --add-stop-codons')
    parser.add_argument('-o', '--out', type=str,
                        help=(f'output GFF3 file, gzipped where the input is (see --gzip) (default: the GFF3 file\'s '
                              f'name with {OUT_SUFFIX}, or {NO_UTRS_SUFFIX} with --no-utrs, before the extension, '
                              f'e.g. genes{OUT_SUFFIX}.gff3)'))
    parser.add_argument('--no-utrs', action='store_true',
                        help=('leave out the exon and UTR lines too, keeping only the CDS lines of each transcript '
                              '(gffcompare then takes the CDS lines as exons)'))
    parser.add_argument('--add-stop-codons', action='store_true',
                        help=('add a stop codon missing from a CDS as GeenuFF\'s import adds it, before picking the '
                              'longest transcript (needs --fasta)'))
    add_output_compression_args(parser)
    parser.add_argument('--log-file', type=str, help='file to write the log to, in addition to the terminal')
    args = parser.parse_args()
    if args.add_stop_codons and not args.fasta:
        parser.error('--add-stop-codons needs --fasta')
    require_input_files(parser, [args.gff, args.fasta])
    require_fasta_and_gff(parser, args.fasta, args.gff)
    default_out = insert_before_extension(args.gff, NO_UTRS_SUFFIX if args.no_utrs else OUT_SUFFIX)
    args.out = output_path(args.out or default_out, output_compression(args, args.gff))
    check_output_paths(parser, [args.out, args.log_file], [args.gff, args.fasta])
    main(args)
