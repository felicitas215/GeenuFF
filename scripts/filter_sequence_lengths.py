#! /usr/bin/env python3
"""Removes short 'junk' sequences from a FASTA file and, given one, their lines from its GFF3 file,
before import or, the FASTA file alone, before prediction: sequences in a statistically distinct
cluster of small ones, e.g. organelle genomes or unplaced scaffolds beside the chromosomes, and,
either way, every sequence below a hard floor (--min-seq-len, by default 1000 bp) (for more see
docs/scripts.md).

The output files default to the input files' names with _length_filtered before the extension,
e.g. genome.fna.gz -> genome_length_filtered.fa.gz. The FASTA file is read in fixed-size chunks,
so a sequence written on a single line is never held in memory whole.
"""
import argparse
import gzip
import logging
import math
import os
from collections import OrderedDict, defaultdict
from typing import BinaryIO, TextIO

from dustdas.fastahelper import text_or_gzip_open

from geenuff.base.helpers import (fasta_header_id, fasta_lines, gff_feature_lines, insert_before_extension,
                                  setup_logging)

logger = logging.getLogger(__name__)

MIN_GAP_RANK = 3.0  # how far the largest gap must stand out from the others, as a robust z-score
# below this many sequences, the gap statistics (median, MAD) are unreliable or degenerate
MIN_SEQUENCES_FOR_CLUSTERING = 5
DEFAULT_MIN_SEQ_LEN = 1000
SEQUENCE_WHITESPACE = b' \t\r\n'
MAX_LOGGED_IDS = 10
OUT_SUFFIX = '_length_filtered'  # inserted before the extension of the input files' names by default


def fasta_lengths(path: str) -> OrderedDict[str, int]:
    """{sequence ID (the header up to the first whitespace): length in bases}, in file order."""
    lengths, seq_id = OrderedDict(), None
    for at_line_start, piece in fasta_lines(path):
        if at_line_start and piece.startswith(b'>'):
            seq_id = fasta_header_id(piece[1:].decode())
            lengths[seq_id] = 0
        elif seq_id is not None:
            lengths[seq_id] += len(piece.translate(None, SEQUENCE_WHITESPACE))
    return lengths


def gap_cluster_ids(lengths: OrderedDict[str, int]) -> set[str]:
    """The IDs of the sequences below the largest gap between the sorted log10 lengths, where that gap
    stands out from the others; empty with too few sequences or no such gap."""
    if len(lengths) < MIN_SEQUENCES_FOR_CLUSTERING:
        return set()
    logs = {seq_id: math.log10(max(length, 1)) for seq_id, length in lengths.items()}
    sorted_logs = sorted(logs.values())
    gaps = [b - a for a, b in zip(sorted_logs, sorted_logs[1:])]
    max_i = max(range(len(gaps)), key=gaps.__getitem__)
    others = sorted(gaps[:max_i] + gaps[max_i + 1:])
    median = others[len(others) // 2]
    mad = sorted(abs(gap - median) for gap in others)[len(others) // 2]
    # where nearly all lengths are identical (very, very unrealistic for a real assembly), the other
    # gaps spread by 0 and any gap, however small, stands out
    if (gaps[max_i] - median) / max(mad, 1e-9) < MIN_GAP_RANK:
        return set()
    threshold = (sorted_logs[max_i] + sorted_logs[max_i + 1]) / 2
    return {seq_id for seq_id, log_length in logs.items() if log_length <= threshold}


def write_fasta_without(fasta_path: str, out_path: str, removed: set[str]) -> None:
    """Copies the FASTA file without the sequences of the removed IDs. A sequence line longer than
    helpers.FASTA_CHUNK_SIZE comes out wrapped into lines of up to that length, which FASTA allows."""
    keep = True
    with open_out(out_path, binary=True) as out:
        for at_line_start, piece in fasta_lines(fasta_path):
            if at_line_start and piece.startswith(b'>'):
                keep = fasta_header_id(piece[1:].decode()) not in removed
            if keep:
                out.write(piece + b'\n')


def write_gff_without(gff_path: str, out_path: str, removed: set[str]) -> None:
    """Copies the GFF3 file without the feature and ##sequence-region lines of the removed IDs."""
    with text_or_gzip_open(gff_path) as f_in, open_out(out_path) as f_out:
        for line in f_in:
            if line.startswith(('##FASTA', '>')):
                break
            if line.startswith('##sequence-region'):
                seqid = line.split()[1]
            elif line.strip() and not line.startswith('#'):
                seqid = line.split('\t', 1)[0]
            else:
                seqid = None
            if seqid not in removed:
                f_out.write(line)


def open_out(path: str, binary: bool = False) -> TextIO | BinaryIO:
    if path.endswith('.gz'):
        return gzip.open(path, 'wb' if binary else 'wt', compresslevel=6)
    return open(path, 'wb' if binary else 'w')


def coding_by_seqid(gff_path: str) -> dict[str, tuple[int, int]]:
    """{sequence ID: (coding genes, coding transcripts) on it}. A coding transcript is a Parent of CDS
    lines and its gene the top of its Parent chain, so a CDS line naming its gene directly counts the
    gene as its own transcript."""
    parent_of, cds_parents = {}, defaultdict(set)
    for line in gff_feature_lines(gff_path):
        columns = line.split('\t')
        if columns[2] == 'exon':  # never a parent
            continue
        attributes = dict(a.split('=', 1) for a in columns[8].strip().split(';') if '=' in a)
        parents = attributes['Parent'].split(',') if 'Parent' in attributes else []
        if columns[2] == 'CDS':
            cds_parents[columns[0]].update(parents)
        elif 'ID' in attributes and parents:
            parent_of[attributes['ID']] = parents[0]

    def gene_of(feature_id: str) -> str:
        seen = set()
        while feature_id in parent_of and feature_id not in seen:
            seen.add(feature_id)
            feature_id = parent_of[feature_id]
        return feature_id

    return {seqid: (len({gene_of(t) for t in ids}), len(ids)) for seqid, ids in cds_parents.items()}


def share(part: int, whole: int) -> str:
    """'part of whole (x.xx%)'"""
    return f'{part} of {whole} ({100 * part / whole if whole else 0:.2f}%)'


def main(args: argparse.Namespace) -> None:
    setup_logging(args.log_file)
    names =' and '.join(f'"{os.path.basename(path)}"' for path in (args.fasta, args.gff) if path)
    logger.info(f'Length filtering {names}')
    lengths = fasta_lengths(args.fasta)
    coding = coding_by_seqid(args.gff) if args.gff else {}
    clustered = set() if args.no_adaptive else gap_cluster_ids(lengths)
    below_floor = {seq_id for seq_id, length in lengths.items() if length < args.min_seq_len}
    to_remove = clustered | below_floor
    removed = [seq_id for seq_id in lengths if seq_id in to_remove]
    logger.info(f'{len(lengths)} sequences: {len(clustered)} in a cluster of small ones'
                + (' (adaptive filter off)' if args.no_adaptive else '')
                + f', {len(below_floor)} below {args.min_seq_len} bp')
    if not removed:
        logger.info('No sequence to remove: nothing written')
        return
    if len(removed) == len(lengths):
        logger.warning('The filter would leave no FASTA sequence, so nothing is removed and nothing written')
        return
    if coding and not set(coding) - to_remove:
        logger.warning('The filter would leave no sequence with coding transcripts, so nothing is removed and '
                       'nothing written')
        return

    removed_lengths = [lengths[seq_id] for seq_id in removed]
    held = ''
    if args.gff:
        lost = [coding[seq_id] for seq_id in removed if seq_id in coding]
        genes, transcripts = sum(g for g, _ in coding.values()), sum(t for _, t in coding.values())
        lost_genes, lost_transcripts = sum(g for g, _ in lost), sum(t for _, t in lost)
        held = (f', and hold {share(lost_genes, genes)} coding genes and {share(lost_transcripts, transcripts)} '
                f'coding transcripts')
    more = len(removed) - MAX_LOGGED_IDS
    logger.info(f'Removing {len(removed)} of {len(lengths)} sequences, {len(lengths) - len(removed)} remaining. '
                f'The removed ones are {min(removed_lengths)}-{max(removed_lengths)} bp long, '
                f'{sum(removed_lengths)} of {sum(lengths.values())} bp{held}; removed: '
                f'{", ".join(removed[:MAX_LOGGED_IDS])}' + (f', ... and {more} more' if more > 0 else ''))
    out_paths = [args.fasta_out or insert_before_extension(args.fasta, OUT_SUFFIX)]
    write_fasta_without(args.fasta, out_paths[0], set(removed))
    if args.gff:
        out_paths.append(args.gff_out or insert_before_extension(args.gff, OUT_SUFFIX))
        write_gff_without(args.gff, out_paths[1], set(removed))
    if args.removed_out:
        with open(args.removed_out, 'w') as f:
            f.writelines(f'{seq_id}\t{lengths[seq_id]}\n' for seq_id in removed)
    written = ' and '.join(f'"{path}"' for path in out_paths)
    logger.info(f'Wrote {written}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--fasta', type=str, required=True, help='FASTA file (may be gzipped)')
    parser.add_argument('--gff', type=str,
                        help='GFF3 file of the FASTA file (may be gzipped); left out to filter the FASTA file alone, '
                             'e.g. before prediction')
    parser.add_argument('--fasta-out', type=str,
                        help=('filtered FASTA file, gzipped if it ends in .gz (default: the FASTA file\'s name with '
                              '_length_filtered before the extension, e.g. genome_length_filtered.fa.gz); not '
                              'written when nothing is removed'))
    parser.add_argument('--gff-out', type=str,
                        help=('filtered GFF3 file, gzipped if it ends in .gz (default: the GFF3 file\'s name with '
                              '_length_filtered before the extension, e.g. genes_length_filtered.gff3); not '
                              'written when nothing is removed'))
    parser.add_argument('--removed-out', type=str, help='TSV file listing each removed sequence ID and its length')
    parser.add_argument('--min-seq-len', type=int, default=DEFAULT_MIN_SEQ_LEN,
                        help=f'hard floor: sequences shorter than this are always removed (default: '
                             f'{DEFAULT_MIN_SEQ_LEN} bp)')
    parser.add_argument('--no-adaptive', action='store_true',
                        help='apply the hard floor alone, without the cluster of small sequences')
    parser.add_argument('--log-file', type=str, help='file to write the log to, in addition to the terminal')
    args = parser.parse_args()
    if args.gff_out and not args.gff:
        parser.error('--gff-out needs --gff')
    main(args)
