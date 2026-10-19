#! /usr/bin/env python3
"""Adds a transcript line for the genes of a GFF3 file that have none, their exon, CDS and UTR lines
naming the gene itself as Parent, as in some MAKER files. GeenuFF masks CDS lines without a
transcript line rather than guess a gene model from them. Each such gene gets an mRNA line spanning
it, written after the gene line with the gene ID and _mRNA as its ID, which the gene's exon, CDS and
UTR lines get as Parent in place of the gene. A gene with transcript lines is left as it is. The log
warns where a gene has overlapping exon or CDS lines, a sign of several isoforms written as one
gene, which GeenuFF masks.
This should only be used on files where each gene parented feature belongs to one isoform only, as
multi-isoform features all parented to one gene can not be sorted into the correct isoform reliably.
"""
import argparse
import logging
import os
from collections import defaultdict

from dustdas.fastahelper import text_or_gzip_open

from geenuff.applications.importer import GENE_TYPES, PIECE_KEYS, TRANSCRIPT_TYPES
from geenuff.base.helpers import (add_output_compression_args, check_output_paths, gff_attributes, gff_feature_lines,
                                  open_out, output_compression, output_path, require_input_files, setup_logging)

logger = logging.getLogger(__name__)

TRANSCRIPT_SUFFIX = '_mRNA'


def overlapping(spans: list[tuple[int, int]]) -> bool:
    spans = sorted(spans)
    return any(start <= previous_end for (_, previous_end), (start, _) in zip(spans, spans[1:]))


def genes_without_transcripts(gff_path: str) -> tuple[set[str], set[str]]:
    """The IDs of the genes without transcript lines whose exon, CDS or UTR lines name them as
    Parent, and of those the ones with overlapping exon or CDS lines.

    Raises:
        SystemExit: where the ID a new transcript line would get is the ID of a line already
    """
    gene_ids, all_ids, with_transcripts = set(), set(), set()
    spans = defaultdict(lambda: defaultdict(list))  # {parent ID: {'exons', 'cds' or 'utrs': [(start, end)]}}
    for line in gff_feature_lines(gff_path):
        columns = line.rstrip('\n').split('\t')
        feature_type, attributes = columns[2], gff_attributes(columns[8])
        all_ids.add(attributes.get('ID'))
        parents = attributes['Parent'].split(',') if 'Parent' in attributes else []
        if feature_type in GENE_TYPES:
            gene_ids.add(attributes.get('ID'))
        elif feature_type in TRANSCRIPT_TYPES:
            with_transcripts.update(parents)
        elif feature_type in PIECE_KEYS:
            for parent in parents:
                spans[parent][PIECE_KEYS[feature_type]].append((int(columns[3]), int(columns[4])))

    genes = (spans.keys() & gene_ids) - with_transcripts
    n_mixed = len((spans.keys() & gene_ids) & with_transcripts)
    if n_mixed:
        logger.info(f'{n_mixed} genes with transcript lines also have exon, CDS or UTR lines naming the gene '
                    f'itself; they are left as they are')
    clashing = sorted(gene_id for gene_id in genes if f'{gene_id}{TRANSCRIPT_SUFFIX}' in all_ids)
    if clashing:
        raise SystemExit(f'{len(clashing)} IDs a new transcript line would get are the IDs of lines already, e.g. '
                         f'{", ".join(gene_id + TRANSCRIPT_SUFFIX for gene_id in clashing[:5])}; no file written')
    overlaps = {gene_id for gene_id in genes if overlapping(spans[gene_id]['exons'])
                or overlapping(spans[gene_id]['cds'])}
    return genes, overlaps


def write_with_transcripts(gff_path: str, out_path: str, genes: set[str], compression_level: int | None) -> None:
    """Writes the GFF3 file with an mRNA line after each gene line of the given IDs, spanning the gene,
    and that mRNA as Parent of the gene's exon, CDS and UTR lines in place of the gene."""
    with text_or_gzip_open(gff_path) as f_in, open_out(out_path, compression_level) as f_out:
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
                if any(parent in genes for parent in parents):
                    new = ','.join(f'{parent}{TRANSCRIPT_SUFFIX}' if parent in genes else parent for parent in parents)
                    columns[8] = columns[8].replace(f'Parent={attributes["Parent"]}', f'Parent={new}')
                    line = '\t'.join(columns) + '\n'
            f_out.write(line)
            gene_id = attributes.get('ID')
            if columns[2] in GENE_TYPES and gene_id in genes:
                seqid, source, _, start, end, _, strand = columns[:7]
                f_out.write(f'{seqid}\t{source}\tmRNA\t{start}\t{end}\t.\t{strand}\t.\t'
                            f'ID={gene_id}{TRANSCRIPT_SUFFIX};Parent={gene_id}\n')


def main(args: argparse.Namespace) -> None:
    setup_logging(args.log_file)
    logger.info(f'Adding transcript lines to the genes without one in "{os.path.basename(args.gff)}"')
    genes, overlaps = genes_without_transcripts(args.gff)
    if not genes:
        logger.info('No gene without transcript lines has exon, CDS or UTR lines naming it: nothing to add, '
                    'no file written')
        return
    logger.info(f'Found {len(genes)} genes without transcript lines whose exon, CDS or UTR lines name the gene, '
                f'e.g. {", ".join(sorted(genes)[:3])}')
    if overlaps:
        logger.warning(f'{len(overlaps)} of them have overlapping exon or CDS lines, a sign of several isoforms '
                       f'written as one gene, which GeenuFF masks, e.g. {", ".join(sorted(overlaps)[:5])}')
    write_with_transcripts(args.gff, args.out, genes, output_compression(args, args.gff))
    logger.info(f'Wrote {len(genes)} transcript lines to "{args.out}"')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--gff', type=str, required=True, help='GFF3 file (may be gzipped)')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help=('output GFF3 file, gzipped where the input is (see --gzip); not written when nothing '
                              'is to be added'))
    add_output_compression_args(parser)
    parser.add_argument('--log-file', type=str, help='file to write the log to, in addition to the terminal')
    args = parser.parse_args()
    require_input_files(parser, [args.gff])
    args.out = output_path(args.out, output_compression(args, args.gff))
    check_output_paths(parser, [args.out, args.log_file], [args.gff])
    main(args)
