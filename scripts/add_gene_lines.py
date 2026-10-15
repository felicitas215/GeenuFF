#! /usr/bin/env python3
"""Adds a gene line for the transcripts of a GFF3 file that have no Parent, as written by gene
predictors such as AUGUSTUS, BRAKER and Tiberius; GeenuFF does not import a transcript without a
gene. A transcript's gene ID is taken from the first of the given attributes it has (by default
geneID, gene_id, gene), else from its ID up to the last occurrence of the split text (by default
".", so g653.t1 gives g653), else it gets a gene of its own, named after it with _gene appended.
The gene line spans the transcripts of its ID, which each get it as their Parent; a gene ID naming
a gene line already is only set as their Parent. The log shows example splits and warns where
transcripts put into one gene do not overlap, a sign of the wrong split text.
"""
import argparse
import logging
import os
from collections import defaultdict
from dataclasses import dataclass, field

from dustdas.fastahelper import text_or_gzip_open

from geenuff.applications.importer import GENE_TYPES, TRANSCRIPT_TYPES
from geenuff.base.helpers import gff_attributes, gff_feature_lines, open_text_out, setup_logging

logger = logging.getLogger(__name__)

DEFAULT_ATTRIBUTES = 'geneID,gene_id,gene'
DEFAULT_SPLIT = '.'


@dataclass
class Gene:
    gene_id: str
    seqid: str
    source: str
    strand: str
    has_line: bool = False  # whether the file has a gene line of this ID already
    transcripts: list[tuple[str, int, int]] = field(default_factory=list)  # (ID, start, end)

    def line(self) -> str:
        start = min(start for _, start, _ in self.transcripts)
        end = max(end for _, _, end in self.transcripts)
        return f'{self.seqid}\t{self.source}\tgene\t{start}\t{end}\t.\t{self.strand}\t.\tID={self.gene_id}\n'

    def overlapping(self) -> bool:
        """Whether its transcripts form one stretch of overlapping transcripts, as isoforms do."""
        spans = sorted((start, end) for _, start, end in self.transcripts)
        reach = spans[0][1]
        for start, end in spans[1:]:
            if start > reach:
                return False
            reach = max(reach, end)
        return True


def gene_id_of(attributes: dict[str, str], gene_attributes: list[str], split_at: str) -> tuple[str, str]:
    """The gene ID of a transcript and where it came from: 'attribute <name>', 'split ID' or 'own ID'."""
    for name in gene_attributes:
        if attributes.get(name):
            return attributes[name], f'attribute {name}'
    transcript_id = attributes['ID']
    cut = transcript_id.rfind(split_at) if split_at else -1
    if cut > 0:
        return transcript_id[:cut], 'split ID'
    return f'{transcript_id}_gene', 'own ID'


def genes_for_orphans(gff_path: str, gene_attributes: list[str], split_at: str) -> dict[str, Gene]:
    """{transcript ID: its gene} for the transcripts without a Parent, their genes named as in the module
    docstring. A gene ID naming a line of another type, or transcripts on several sequences or strands,
    gets the transcripts' sequence and strand appended to stay unique."""
    gene_line_ids, other_ids = set(), set()
    orphans = []  # (transcript ID, gene ID, seqid, source, strand, start, end)
    sources, split_examples = defaultdict(int), []
    n_without_id = 0  # transcripts with neither Parent nor ID, which no line can name as Parent
    for line in gff_feature_lines(gff_path):
        seqid, source, feature_type, start, end, _, strand, _, column = line.rstrip('\n').split('\t')[:9]
        attributes = gff_attributes(column)
        line_id = attributes.get('ID')
        (gene_line_ids if feature_type in GENE_TYPES else other_ids).add(line_id)
        if feature_type in TRANSCRIPT_TYPES and 'Parent' not in attributes and not line_id:
            n_without_id += 1
        elif feature_type in TRANSCRIPT_TYPES and 'Parent' not in attributes:
            gene_id, how = gene_id_of(attributes, gene_attributes, split_at)
            sources[how] += 1
            if how == 'split ID' and len(split_examples) < 3:
                split_examples.append(f'{line_id} -> {gene_id}')
            orphans.append((line_id, gene_id, seqid, source, strand, int(start), int(end)))

    places = defaultdict(set)
    for _, gene_id, seqid, _, strand, _, _ in orphans:
        places[gene_id].add((seqid, strand))
    several_places = {gene_id for gene_id, where in places.items() if len(where) > 1}
    clashing = {gene_id for gene_id in places if gene_id in other_ids and gene_id not in gene_line_ids}
    genes, by_transcript = {}, {}
    for transcript_id, gene_id, seqid, source, strand, start, end in orphans:
        if gene_id in several_places or gene_id in clashing:
            gene_id = f'{gene_id}_{seqid}_{"plus" if strand == "+" else "minus"}'
        gene = genes.setdefault(gene_id, Gene(gene_id, seqid, source, strand, has_line=gene_id in gene_line_ids))
        gene.transcripts.append((transcript_id, start, end))
        by_transcript[transcript_id] = gene

    if n_without_id:
        logger.warning(f'{n_without_id} transcripts have neither a Parent nor an ID and are left as they are; no '
                       f'exon or CDS line can name them as Parent')
    if orphans:
        logger.info(f'Found {len(orphans)} transcripts without a Parent, their gene IDs taken from '
                    + ', '.join(f'{how} for {n}' for how, n in sorted(sources.items(), key=lambda item: -item[1])))
    if split_examples:
        logger.info(f'  split at the last "{split_at}", e.g. {", ".join(split_examples)}')
    if sources['own ID']:
        logger.warning(f'{sources["own ID"]} transcripts have none of the attributes {", ".join(gene_attributes)} '
                       f'nor "{split_at}" in their ID, so each gets a gene of its own; set --gene-id-attributes or '
                       f'--split-at if isoforms should share one')
    apart = sorted(gene.gene_id for gene in genes.values() if not gene.overlapping())
    if apart:
        logger.warning(f'{len(apart)} genes hold transcripts that do not overlap each other, which isoforms do; '
                       f'check --split-at and --gene-id-attributes, e.g. {", ".join(apart[:5])}')
    if several_places:
        logger.warning(f'{len(several_places)} gene IDs name transcripts on several sequences or strands, split into '
                       f'one gene per sequence and strand, e.g. {", ".join(sorted(several_places)[:5])}')
    if clashing:
        logger.warning(f'{len(clashing)} gene IDs are the ID of another line already, given the sequence and strand '
                       f'as a suffix, e.g. {", ".join(sorted(clashing)[:5])}')
    n_existing = sum(gene.has_line for gene in genes.values())
    if n_existing:
        logger.info(f'{n_existing} gene IDs name a gene line already, which their transcripts get as Parent')
    return by_transcript


def write_with_genes(gff_path: str, out_path: str, gene_of: dict[str, Gene]) -> int:
    """Writes the GFF3 file with each new gene line before its first transcript and each orphan
    transcript's gene as its Parent. Returns the number of gene lines written."""
    written = set()
    with text_or_gzip_open(gff_path) as f_in, open_text_out(out_path) as f_out:
        for line in f_in:
            if line.startswith(('##FASTA', '>')):
                break
            columns = line.rstrip('\n').split('\t')
            if len(columns) >= 9 and columns[2] in TRANSCRIPT_TYPES:
                attributes = gff_attributes(columns[8])
                gene = gene_of.get(attributes.get('ID'))
                if gene is not None and 'Parent' not in attributes:
                    if not gene.has_line and gene.gene_id not in written:
                        f_out.write(gene.line())
                        written.add(gene.gene_id)
                    columns[8] = f'{columns[8].strip().rstrip(";")};Parent={gene.gene_id}'
                    line = '\t'.join(columns) + '\n'
            f_out.write(line)
    return len(written)


def main(args: argparse.Namespace) -> None:
    setup_logging(args.log_file)
    logger.info(f'Adding gene lines to the transcripts without a Parent in "{os.path.basename(args.gff)}"')
    gene_attributes = [name for name in args.gene_id_attributes.split(',') if name]
    gene_of = genes_for_orphans(args.gff, gene_attributes, args.split_at)
    if not gene_of:
        logger.info('Every transcript with an ID has a Parent: nothing to add, no file written')
        return
    n_genes = write_with_genes(args.gff, args.out, gene_of)
    logger.info(f'Wrote {n_genes} gene lines for {len(gene_of)} transcripts to "{args.out}"')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--gff', type=str, required=True, help='GFF3 file (may be gzipped)')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output GFF3 file, gzipped if it ends in .gz; not written when nothing is to be added')
    parser.add_argument('--gene-id-attributes', type=str, default=DEFAULT_ATTRIBUTES,
                        help=(f'comma-separated attributes of a transcript line holding its gene ID, the first one '
                              f'present used (default: {DEFAULT_ATTRIBUTES}); give an empty value, written '
                              f'--gene-id-attributes "" or --gene-id-attributes=, to split the transcript ID alone'))
    parser.add_argument('--split-at', type=str, default=DEFAULT_SPLIT,
                        help=(f'text at whose last occurrence a transcript ID is cut to get its gene ID, where no '
                              f'gene ID attribute is present (default: "{DEFAULT_SPLIT}"). E.g. "." turns g653.t1 '
                              f'into g653, and "-mRNA-" turns gene7-mRNA-2 into gene7'))
    parser.add_argument('--log-file', type=str, help='file to write the log to, in addition to the terminal')
    main(parser.parse_args())
