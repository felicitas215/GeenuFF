#! /usr/bin/env python3
"""Removes short 'junk' sequences from a FASTA file and, given one, their lines from its GFF3 file,
before import or, the FASTA file alone, before prediction: sequences in a statistically distinct
cluster of small ones, e.g. organelle genomes or unplaced scaffolds beside the chromosomes, and,
either way, every sequence below a hard floor (--min-seq-len, by default 1000 bp) (for more see
docs/scripts.md).

The output files default to the input files' names with _length_filtered before the extension,
e.g. genome.fna.gz -> genome_length_filtered.fna.gz, each gzipped where its input is (see --gzip).
The FASTA file is read in fixed-size chunks, so a sequence written on a single line is never held
in memory whole.
"""
import argparse
import logging
import math
import os
from collections import OrderedDict, defaultdict

from geenuff.base.helpers import (add_output_compression_args, check_output_paths, insert_before_extension,
                                  output_compression, output_path, require_fasta_and_gff, require_input_files,
                                  setup_logging, write_fasta_without, write_gff_without)
from geenuff.base.types import PIECE_KEYS
from geenuff.parsing.fastahelper import SEQUENCE_WHITESPACE, FastaParser
from geenuff.parsing.gffhelper import GFFObject, GFFParser

logger = logging.getLogger(__name__)

MIN_GAP_RANK = 3.0  # how far the largest gap must stand out from the others, as a robust z-score
# how many times shorter the longest sequence of the cluster must be than the shortest one kept. In a
# chromosome-level assembly all gaps are small, so one barely larger stands out and would cut off the
# smallest chromosome; the cuts at real chromosomes seen were 1.19- to 2.42-fold steps, those at
# organelles and junk scaffolds 2.68-fold and more
MIN_CLUSTER_STEP = 2.5
# below this many sequences, the gap statistics (median, MAD) are unreliable or degenerate
MIN_SEQUENCES_FOR_CLUSTERING = 5
# junk sequences hold a small share of an assembly's length; a larger cluster means the gap lies within one
# continuous length distribution, as in fragmented scaffold-level assemblies, and nothing is removed
DEFAULT_MAX_CLUSTER_LENGTH_FRACTION = 0.05
DEFAULT_MIN_SEQ_LEN = 1000
# above this share of the bases or coding genes, sequences below the floor are still removed, with a warning
FLOOR_WARN_SHARE = 0.05
MAX_LOGGED_IDS = 10
OUT_SUFFIX = '_length_filtered'  # inserted before the extension of the input files' names by default


def fasta_lengths(path: str) -> OrderedDict[str, int]:
    """{sequence ID (the header up to the first whitespace): length in bases}, in file order."""
    lengths = OrderedDict()
    for is_header, piece in FastaParser.pieces(path):
        if is_header:
            seq_id = FastaParser.header_id(piece[1:].decode())
            lengths[seq_id] = 0
        else:
            lengths[seq_id] += len(piece.translate(None, SEQUENCE_WHITESPACE))
    return lengths


def gap_cluster_ids(lengths: OrderedDict[str, int], max_length_fraction: float) -> set[str]:
    """The IDs of the sequences below the largest gap between the sorted log10 lengths, where that gap
    stands out from the others and is a step of at least MIN_CLUSTER_STEP in length; empty with too few
    sequences, no such gap or a cluster holding more than max_length_fraction of the total length."""
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
    step = 10 ** gaps[max_i]
    if step < MIN_CLUSTER_STEP:
        logger.info(f'The largest gap between the sequence lengths, at {round(10 ** sorted_logs[max_i])} bp, stands '
                    f'out but is only a {step:.2f}-fold step, under the {MIN_CLUSTER_STEP:g}-fold one separating '
                    f'junk from real sequences, e.g. the smallest chromosome, so no cluster is removed')
        return set()
    threshold = (sorted_logs[max_i] + sorted_logs[max_i + 1]) / 2
    cluster = {seq_id for seq_id, log_length in logs.items() if log_length <= threshold}
    cluster_length, total_length = sum(lengths[seq_id] for seq_id in cluster), sum(lengths.values())
    if cluster_length > max_length_fraction * total_length:
        logger.warning(f'The cluster of {len(cluster)} small sequences holds '
                       f'{100 * cluster_length / total_length:.2f}% of the total length ({cluster_length} of '
                       f'{total_length} bp), more than the {100 * max_length_fraction:g}% allowed '
                       f'(--max-cluster-length-fraction {max_length_fraction}), so it is kept')
        return set()
    return cluster


def coding_by_seqid(gff_path: str) -> dict[str, tuple[int, int]]:
    """{sequence ID: (coding genes, coding transcripts) on it}. A coding transcript is a Parent of CDS
    lines and its gene the top of its Parent chain, so a CDS line naming its gene directly counts the
    gene as its own transcript."""
    parent_of, cds_parents = {}, defaultdict(set)
    for line in GFFParser.feature_lines(gff_path):
        entry = GFFObject(line)
        parents = entry.get_Parent() or []
        piece = PIECE_KEYS.get(entry.type)
        if piece == 'cds':
            cds_parents[entry.seqid].update(parents)
        # exon and UTR lines are never a parent, so only gene and transcript lines make the chain
        elif piece is None and entry.get_ID() is not None and parents:
            parent_of[entry.get_ID()] = parents[0]

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
    names = ' and '.join(f'"{os.path.basename(path)}"' for path in (args.fasta, args.gff) if path)
    logger.info(f'Length filtering {names}')
    lengths = fasta_lengths(args.fasta)
    coding = coding_by_seqid(args.gff) if args.gff else {}
    clustered = set() if args.no_adaptive else gap_cluster_ids(lengths, args.max_cluster_length_fraction)
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

    floor_bp_share = sum(lengths[seq_id] for seq_id in below_floor) / sum(lengths.values())
    all_genes = sum(g for g, _ in coding.values())
    floor_gene_share = (sum(coding[seq_id][0] for seq_id in below_floor if seq_id in coding) / all_genes
                        if all_genes else 0)
    if max(floor_bp_share, floor_gene_share) > FLOOR_WARN_SHARE:
        logger.warning(f'The {len(below_floor)} sequences below {args.min_seq_len} bp hold '
                       f'{100 * floor_bp_share:.2f}% of the bases'
                       + (f' and {100 * floor_gene_share:.2f}% of the coding genes' if all_genes else '')
                       + f', more than {100 * FLOOR_WARN_SHARE:g}%. They are removed, as sequences this short '
                         f'rarely give good training data or predictions; a lower --min-seq-len (0 for none) '
                         f'keeps them')

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
    out_paths = [args.fasta_out]
    write_fasta_without(args.fasta, args.fasta_out, set(removed), output_compression(args, args.fasta))
    if args.gff:
        out_paths.append(args.gff_out)
        write_gff_without(args.gff, args.gff_out, set(removed), output_compression(args, args.gff))
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
                        help=('filtered FASTA file, gzipped where the FASTA file is (see --gzip) (default: the FASTA '
                              'file\'s name with _length_filtered before the extension, e.g. '
                              'genome_length_filtered.fa.gz); not written when nothing is removed'))
    parser.add_argument('--gff-out', type=str,
                        help=('filtered GFF3 file, gzipped where the GFF3 file is (see --gzip) (default: the GFF3 '
                              'file\'s name with _length_filtered before the extension, e.g. '
                              'genes_length_filtered.gff3); not written when nothing is removed'))
    parser.add_argument('--removed-out', type=str, help='TSV file listing each removed sequence ID and its length')
    parser.add_argument('--min-seq-len', type=int, default=DEFAULT_MIN_SEQ_LEN,
                        help=f'hard floor: sequences shorter than this are always removed (default: '
                             f'{DEFAULT_MIN_SEQ_LEN} bp)')
    parser.add_argument('--no-adaptive', action='store_true',
                        help='apply the hard floor alone, without the cluster of small sequences')
    parser.add_argument('--max-cluster-length-fraction', type=float, default=DEFAULT_MAX_CLUSTER_LENGTH_FRACTION,
                        help=f'the cluster of small sequences is kept where it holds more than this share of the '
                             f'total length (default: {DEFAULT_MAX_CLUSTER_LENGTH_FRACTION})')
    add_output_compression_args(parser)
    parser.add_argument('--log-file', type=str, help='file to write the log to, in addition to the terminal')
    args = parser.parse_args()
    require_input_files(parser, [args.fasta, args.gff])
    require_fasta_and_gff(parser, args.fasta, args.gff)
    if args.gff_out and not args.gff:
        parser.error('--gff-out needs --gff')
    if not 0 <= args.max_cluster_length_fraction <= 1:
        parser.error('--max-cluster-length-fraction must be between 0 and 1')
    args.fasta_out = output_path(args.fasta_out or insert_before_extension(args.fasta, OUT_SUFFIX),
                                 output_compression(args, args.fasta))
    if args.gff:
        args.gff_out = output_path(args.gff_out or insert_before_extension(args.gff, OUT_SUFFIX),
                                   output_compression(args, args.gff))
    check_output_paths(parser, [args.fasta_out, args.gff_out, args.removed_out, args.log_file],
                       [args.fasta, args.gff])
    main(args)
