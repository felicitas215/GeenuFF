#! /usr/bin/env python3
"""Rebuilds the exon lines of the transcripts of a GFF3 file whose CDS lines reach outside them, as
files can rarely list exon lines for the UTRs only. The new exons are the union of the transcript's exon,
CDS and UTR lines, touching or overlapping lines forming one exon; they replace its old exon lines and
follow its transcript line. Other transcripts are left as they are.

Only for files where this is systematic: elsewhere, a CDS outside its exons means the CDS or the
exons are wrong, which GeenuFF masks (cds_outside_exons) and rebuilding would hide. The log gives the
share of coding transcripts rebuilt; check a few rebuilt genes against the source.
"""
import argparse
import logging
import os
from collections import defaultdict

from dustdas.fastahelper import text_or_gzip_open

from geenuff.applications.importer import CDS_OUTSIDE_EXONS_WARN_SHARE, PIECE_KEYS, TRANSCRIPT_TYPES
from geenuff.base.helpers import gff_attributes, gff_feature_lines, open_text_out, setup_logging

logger = logging.getLogger(__name__)


def merged(spans: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Sorted GFF spans (1-based, inclusive) merged where they touch or overlap."""
    out = []
    for start, end in sorted(spans):
        if out and start <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], end)
        else:
            out.append([start, end])
    return [(start, end) for start, end in out]


def rebuilt_exons(gff_path: str) -> tuple[dict[str, list[tuple[int, int]]], int]:
    """{transcript ID: its new exon spans} for the transcripts with exon lines whose CDS lines reach
    outside them, and the number of coding transcripts with exon lines."""
    transcript_ids = set()
    pieces = defaultdict(lambda: defaultdict(list))  # {transcript ID: {'exons'/'cds'/'utrs': [spans]}}
    for line in gff_feature_lines(gff_path):
        columns = line.rstrip('\n').split('\t')
        feature_type, attributes = columns[2], gff_attributes(columns[8])
        if feature_type in TRANSCRIPT_TYPES and 'ID' in attributes:
            transcript_ids.add(attributes['ID'])
        elif feature_type in PIECE_KEYS:
            for parent in attributes.get('Parent', '').split(','):
                pieces[parent][PIECE_KEYS[feature_type]].append((int(columns[3]), int(columns[4])))

    rebuilt, n_coding = {}, 0
    for transcript_id, spans in pieces.items():
        if transcript_id not in transcript_ids or not spans['cds'] or not spans['exons']:
            continue
        n_coding += 1
        exons = merged(spans['exons'])
        if any(not any(lo <= start and end <= hi for lo, hi in exons) for start, end in spans['cds']):
            rebuilt[transcript_id] = merged(spans['exons'] + spans['cds'] + spans['utrs'])
    return rebuilt, n_coding


def write_rebuilt(gff_path: str, out_path: str, rebuilt: dict[str, list[tuple[int, int]]]) -> None:
    """Writes the GFF3 file with the old exon lines of the rebuilt transcripts left out (an exon line
    shared with other transcripts keeps those as Parent) and their new exon lines after each
    transcript line."""
    with text_or_gzip_open(gff_path) as f_in, open_text_out(out_path) as f_out:
        for line in f_in:
            if line.startswith(('##FASTA', '>')):
                break
            columns = line.rstrip('\n').split('\t')
            if len(columns) < 9 or line.startswith('#'):
                f_out.write(line)
                continue
            attributes = gff_attributes(columns[8])
            if PIECE_KEYS.get(columns[2]) == 'exons':
                parents = attributes.get('Parent', '').split(',')
                kept = [parent for parent in parents if parent not in rebuilt]
                if not kept:
                    continue
                if len(kept) < len(parents):
                    columns[8] = columns[8].replace(f'Parent={attributes["Parent"]}', f'Parent={",".join(kept)}')
                    line = '\t'.join(columns) + '\n'
            f_out.write(line)
            transcript_id = attributes.get('ID')
            if columns[2] in TRANSCRIPT_TYPES and transcript_id in rebuilt:
                seqid, source, strand = columns[0], columns[1], columns[6]
                for start, end in rebuilt[transcript_id]:
                    f_out.write(f'{seqid}\t{source}\texon\t{start}\t{end}\t.\t{strand}\t.\tParent={transcript_id}\n')


def main(args: argparse.Namespace) -> None:
    setup_logging(args.log_file)
    logger.info(f'Rebuilding the exon lines of transcripts with CDS lines outside them in '
                f'"{os.path.basename(args.gff)}"')
    rebuilt, n_coding = rebuilt_exons(args.gff)
    if not rebuilt:
        logger.info('Every CDS line lies within its exon lines: nothing to rebuild, no file written')
        return
    share = len(rebuilt) / n_coding
    examples = ', '.join(sorted(rebuilt)[:5])
    logger.info(f'Rebuilt the exons of {len(rebuilt)} of {n_coding} coding transcripts with exon lines '
                f'({share:.1%}), e.g. {examples}')
    if share <= CDS_OUTSIDE_EXONS_WARN_SHARE:
        logger.warning(f'Only {share:.1%} of the coding transcripts have CDS lines outside their exons, which '
                       f'looks like a few wrong genes rather than exon lines for the UTRs only; GeenuFF masks '
                       f'those on its own, so the original file is likely the better choice')
    write_rebuilt(args.gff, args.out, rebuilt)
    logger.info(f'Wrote "{args.out}"')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--gff', type=str, required=True, help='GFF3 file (may be gzipped)')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output GFF3 file, gzipped if it ends in .gz; not written when nothing is to be rebuilt')
    parser.add_argument('--log-file', type=str, help='file to write the log to, in addition to the terminal')
    args = parser.parse_args()
    if os.path.realpath(args.out) == os.path.realpath(args.gff):
        parser.error('-o must differ from --gff, which is read while the output is written')
    main(args)
