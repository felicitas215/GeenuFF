#! /usr/bin/env python3
"""Puts the exon, CDS and UTR lines of a GFF3 file that name a gene itself as Parent under a
transcript of that gene, as GeenuFF leaves out or masks them otherwise. A gene without transcript
lines, as in some MAKER files, gets an mRNA line spanning it, written after the gene line with the
gene ID and _mRNA as its ID. A gene whose only transcript has no lines of its own gets them under
that transcript. Either way the transcript becomes their Parent in place of the gene. Genes with
several transcripts, or one with lines of its own, are left as they are. The log warns where a gene
has overlapping exon or CDS lines, a sign of several isoforms written as one gene, which GeenuFF
masks.
This should in general only be used on files where each gene parented feature belongs to one isoform
only, as multi-isoform features all parented to one gene can not be sorted into the correct isoform
reliably; which is why genes with several transcripts are left as they are.
"""
import argparse
import logging
import os
from collections import defaultdict

from dustdas.fastahelper import text_or_gzip_open

from geenuff.applications.importer import GENE_TYPES, PIECE_KEYS, TRANSCRIPT_TYPES
from geenuff.base.helpers import (Compression, add_output_compression_args, check_output_paths, gff_attributes,
                                  gff_feature_lines, insert_before_extension, open_out, output_compression, output_path,
                                  require_fasta_and_gff, require_input_files, setup_logging)

logger = logging.getLogger(__name__)

TRANSCRIPT_SUFFIX = '_mRNA'
OUT_SUFFIX = '_with_transcripts'  # inserted before the extension of the GFF3 file's name for the default output


def overlapping(spans: list[tuple[int, int]]) -> bool:
    spans = sorted(spans)
    return any(start <= previous_end for (_, previous_end), (start, _) in zip(spans, spans[1:]))


def transcripts_for_genes(gff_path: str) -> tuple[dict[str, str], set[str], set[str]]:
    """{gene ID: the transcript its exon, CDS and UTR lines naming the gene go under} for the genes
    without transcript lines and those whose only transcript has no lines of its own; of those, the
    genes getting a new transcript line, and the ones with overlapping exon or CDS lines.

    Raises:
        SystemExit: where the ID a new transcript line would get is the ID of a line already
    """
    gene_ids, all_ids = set(), set()
    transcripts_of = defaultdict(list)  # {gene ID: [its transcript IDs]}
    spans = defaultdict(lambda: defaultdict(list))  # {parent ID: {'exons', 'cds' or 'utrs': [(start, end)]}}
    for line in gff_feature_lines(gff_path):
        columns = line.rstrip('\n').split('\t')
        feature_type, attributes = columns[2], gff_attributes(columns[8])
        all_ids.add(attributes.get('ID'))
        parents = attributes['Parent'].split(',') if 'Parent' in attributes else []
        if feature_type in GENE_TYPES:
            gene_ids.add(attributes.get('ID'))
        elif feature_type in TRANSCRIPT_TYPES:
            for parent in parents:
                transcripts_of[parent].append(attributes.get('ID'))
        elif feature_type in PIECE_KEYS:
            for parent in parents:
                spans[parent][PIECE_KEYS[feature_type]].append((int(columns[3]), int(columns[4])))

    gene_parented = spans.keys() & gene_ids
    new = gene_parented - transcripts_of.keys()
    target = {gene_id: f'{gene_id}{TRANSCRIPT_SUFFIX}' for gene_id in new}
    target |= {gene_id: transcripts_of[gene_id][0] for gene_id in gene_parented - new
               if len(transcripts_of[gene_id]) == 1 and transcripts_of[gene_id][0] not in spans}
    n_left = len(gene_parented) - len(target)
    if n_left:
        logger.info(f'{n_left} genes with several transcripts, or one with lines of its own, also have exon, CDS or '
                    f'UTR lines naming the gene itself; they are left as they are')
    clashing = sorted(target[gene_id] for gene_id in new if target[gene_id] in all_ids)
    if clashing:
        raise SystemExit(f'{len(clashing)} IDs a new transcript line would get are the IDs of lines already, e.g. '
                         f'{", ".join(clashing[:5])}; no file written')
    overlaps = {gene_id for gene_id in target if overlapping(spans[gene_id]['exons'])
                or overlapping(spans[gene_id]['cds'])}
    return target, new, overlaps


def write_with_transcripts(gff_path: str, out_path: str, target: dict[str, str], new: set[str],
                           compression: Compression | None) -> None:
    """Writes the GFF3 file with each gene's target transcript (see transcripts_for_genes) as Parent of
    its exon, CDS and UTR lines in place of the gene, and an mRNA line spanning the gene after each
    gene line of the new ones."""
    with text_or_gzip_open(gff_path) as f_in, open_out(out_path, compression) as f_out:
        for line in f_in:
            if line.startswith(('##FASTA', '>')):
                break
            columns = line.rstrip('\n').split('\t')
            if len(columns) < 9 or line.startswith('#'):
                f_out.write(line)
                continue
            attributes = gff_attributes(columns[8])
            if columns[2] in PIECE_KEYS and 'Parent' in attributes:
                parents = attributes['Parent'].split(',')
                if any(parent in target for parent in parents):
                    renamed = ','.join(target.get(parent, parent) for parent in parents)
                    columns[8] = columns[8].replace(f'Parent={attributes["Parent"]}', f'Parent={renamed}')
                    line = '\t'.join(columns) + '\n'
            f_out.write(line)
            gene_id = attributes.get('ID')
            if columns[2] in GENE_TYPES and gene_id in new:
                seqid, source, _, start, end, _, strand = columns[:7]
                f_out.write(f'{seqid}\t{source}\tmRNA\t{start}\t{end}\t.\t{strand}\t.\t'
                            f'ID={target[gene_id]};Parent={gene_id}\n')


def main(args: argparse.Namespace) -> None:
    setup_logging(args.log_file)
    logger.info(f'Putting the exon, CDS and UTR lines naming a gene under a transcript in '
                f'"{os.path.basename(args.gff)}"')
    target, new, overlaps = transcripts_for_genes(args.gff)
    if not target:
        logger.info('No exon, CDS or UTR line names a gene without transcripts, or one whose only transcript has '
                    'no lines of its own: nothing to change, no file written')
        return
    for genes, which in ((sorted(new), 'without transcript lines, each given an mRNA line'),
                         (sorted(target.keys() - new), 'whose only transcript has no lines of its own')):
        if genes:
            logger.info(f'Found {len(genes)} genes {which}, e.g. {", ".join(genes[:3])}')
    if overlaps:
        logger.warning(f'{len(overlaps)} of them have overlapping exon or CDS lines, a sign of several isoforms '
                       f'written as one gene, which GeenuFF masks, e.g. {", ".join(sorted(overlaps)[:5])}')
    write_with_transcripts(args.gff, args.out, target, new, output_compression(args, args.gff))
    logger.info(f'Wrote "{args.out}": the lines naming {len(target)} genes put under a transcript, '
                f'{len(new)} transcript lines added')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--gff', type=str, required=True, help='GFF3 file (may be gzipped)')
    parser.add_argument('-o', '--out', type=str,
                        help=(f'output GFF3 file, gzipped where the input is (see --gzip) (default: the GFF3 file\'s '
                              f'name with {OUT_SUFFIX} before the extension, e.g. genes{OUT_SUFFIX}.gff3); not '
                              f'written when nothing is to be added'))
    add_output_compression_args(parser)
    parser.add_argument('--log-file', type=str, help='file to write the log to, in addition to the terminal')
    args = parser.parse_args()
    require_input_files(parser, [args.gff])
    require_fasta_and_gff(parser, None, args.gff)
    default_out = insert_before_extension(args.gff, OUT_SUFFIX)
    args.out = output_path(args.out or default_out, output_compression(args, args.gff))
    check_output_paths(parser, [args.out, args.log_file], [args.gff])
    main(args)
