#! /usr/bin/env python3
"""Removes short 'junk' sequences from a FASTA file and their lines from its GFF3 file before import:
sequences in a statistically distinct cluster of small ones, e.g. organelle genomes or unplaced
scaffolds beside the chromosomes, and, either way, every sequence below a hard floor
(--min-seq-len, by default 1000 bp) (for more see docs/scripts.md).

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
from collections.abc import Iterator
from typing import BinaryIO, TextIO

from dustdas.fastahelper import text_or_gzip_open

from geenuff.base.helpers import gff_feature_lines

logger = logging.getLogger(__name__)

MIN_GAP_RANK = 3.0  # how far the largest gap must stand out from the others, as a robust z-score
# below this many sequences, the gap statistics (median, MAD) are unreliable or degenerate
MIN_SEQUENCES_FOR_CLUSTERING = 5
DEFAULT_MIN_SEQ_LEN = 1000
CHUNK_SIZE = 1 << 20
SEQUENCE_WHITESPACE = b' \t\r\n'
MAX_LOGGED_IDS = 10


def open_binary(path: str) -> BinaryIO:
    with open(path, 'rb') as f:
        is_gzip = f.read(2) == b'\x1f\x8b'
    return gzip.open(path, 'rb') if is_gzip else open(path, 'rb')


def fasta_lines(path: str) -> Iterator[tuple[bool, bytes]]:
    """Yields (whether it begins a line, piece) for the pieces of each FASTA line, a line longer than
    CHUNK_SIZE coming in several pieces, without their line breaks."""
    at_line_start, rest = True, b''
    with open_binary(path) as f:
        while chunk := f.read(CHUNK_SIZE):
            lines = (rest + chunk).split(b'\n')
            rest = lines.pop()
            for line in lines:
                yield at_line_start, line
                at_line_start = True
            # a line running on past the chunk is passed on in pieces, its header bytes kept whole
            if rest and not (at_line_start and rest.startswith(b'>')) and len(rest) > CHUNK_SIZE:
                yield at_line_start, rest
                at_line_start, rest = False, b''
    if rest:
        yield at_line_start, rest


def fasta_lengths(path: str) -> OrderedDict[str, int]:
    """{sequence ID (the header up to the first whitespace): length in bases}, in file order."""
    lengths, seq_id = OrderedDict(), None
    for at_line_start, piece in fasta_lines(path):
        if at_line_start and piece.startswith(b'>'):
            seq_id = piece[1:].split()[0].decode()
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
    if (gaps[max_i] - median) / max(mad, 1e-9) < MIN_GAP_RANK:
        return set()
    threshold = (sorted_logs[max_i] + sorted_logs[max_i + 1]) / 2
    return {seq_id for seq_id, log_length in logs.items() if log_length <= threshold}


def write_fasta_without(fasta_path: str, out_path: str, removed: set[str]) -> None:
    """Copies the FASTA file without the sequences of the removed IDs. A sequence line longer than
    CHUNK_SIZE comes out wrapped into lines of up to that length, which FASTA allows."""
    keep = True
    with open_out(out_path, binary=True) as out:
        for at_line_start, piece in fasta_lines(fasta_path):
            if at_line_start and piece.startswith(b'>'):
                keep = piece[1:].split()[0].decode() not in removed
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


def default_out_path(path: str) -> str:
    """The path with _length_filtered before its extension, e.g. genome.fna.gz -> genome_length_filtered.fna.gz"""
    gz = '.gz' if path.endswith('.gz') else ''
    stem, extension = os.path.splitext(path[:len(path) - len(gz)])
    return f'{stem}_length_filtered{extension}{gz}'


def open_out(path: str, binary: bool = False) -> TextIO | BinaryIO:
    if path.endswith('.gz'):
        return gzip.open(path, 'wb' if binary else 'wt', compresslevel=6)
    return open(path, 'wb' if binary else 'w')


def coding_transcripts_by_seqid(gff_path: str) -> dict[str, int]:
    """{sequence ID: how many transcripts with CDS lines lie on it}"""
    parents = defaultdict(set)
    for line in gff_feature_lines(gff_path):
        columns = line.split('\t')
        if columns[2] == 'CDS':
            for attribute in columns[8].strip().split(';'):
                if attribute.startswith('Parent='):
                    parents[columns[0]].update(attribute[len('Parent='):].split(','))
    return {seqid: len(ids) for seqid, ids in parents.items()}


def main(args: argparse.Namespace) -> None:
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')
    logger.info(f'Length filtering "{os.path.basename(args.fasta)}" and "{os.path.basename(args.gff)}"')
    lengths = fasta_lengths(args.fasta)
    coding = coding_transcripts_by_seqid(args.gff)
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
    n_coding, n_lost = sum(coding.values()), sum(coding.get(seq_id, 0) for seq_id in removed)
    more = len(removed) - MAX_LOGGED_IDS
    logger.info(f'Removing {len(removed)} of {len(lengths)} sequences, {len(lengths) - len(removed)} remaining. '
                f'The removed ones are {min(removed_lengths)}-{max(removed_lengths)} bp long, '
                f'{sum(removed_lengths)} of {sum(lengths.values())} bp, and hold {n_lost} of {n_coding} '
                f'coding transcripts ({100 * n_lost / n_coding if n_coding else 0:.2f}%); removed: '
                f'{", ".join(removed[:MAX_LOGGED_IDS])}' + (f', ... and {more} more' if more > 0 else ''))
    fasta_out = args.fasta_out or default_out_path(args.fasta)
    gff_out = args.gff_out or default_out_path(args.gff)
    write_fasta_without(args.fasta, fasta_out, set(removed))
    write_gff_without(args.gff, gff_out, set(removed))
    if args.removed_out:
        with open(args.removed_out, 'w') as f:
            f.writelines(f'{seq_id}\t{lengths[seq_id]}\n' for seq_id in removed)
    logger.info(f'Wrote "{fasta_out}" and "{gff_out}"')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--fasta', type=str, required=True, help='FASTA file (may be gzipped)')
    parser.add_argument('--gff', type=str, required=True, help='GFF3 file of the FASTA file (may be gzipped)')
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
    main(parser.parse_args())
