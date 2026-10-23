#! /usr/bin/env python3
"""Lists the transcripts of a GeenuFF database, or their gene or transcript IDs, filtered by error
types, error severity, masked in full or not, exported or not and why not. By default, every coding
transcript is listed, every isoform; filters combine with AND.

Table columns: gene_id, transcript_id, seqid, strand, start, end (GFF3, 1-based), coding,
selected_for_export (the one transcript per gene an export takes), exported (selected and its gene
not excluded), excluded_from_export (why not), errors, error_severity, masked_in_full and masks
(GFF3 ranges). Non-coding transcripts and those of genes left out before the error checks have
"not checked" as errors and match no error filter.
"""
import argparse
import logging
import sqlite3
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field

from geenuff.applications.import_summary import UNEXPORTED_REASONS
from geenuff.applications.importer import error_severity
from geenuff.base import types
from geenuff.base.helpers import check_output_paths, geenuff_to_gff_start_end, is_masked_whole, require_input_files

logger = logging.getLogger(__name__)

ERROR_TYPES = sorted(e.value for e in types.Errors)
SEVERITIES = ['none', 'flank_masked', 'masked_outright']
COLUMNS = ['gene_id', 'transcript_id', 'seqid', 'strand', 'start', 'end', 'coding', 'selected_for_export',
           'exported', 'excluded_from_export', 'errors', 'error_severity', 'masked_in_full', 'masks']
# genes left out for these reasons are excluded before the error checks
UNCHECKED_REASONS = set(types.unrepresentable_reasons)


@dataclass
class Feature:
    start: int
    end: int
    is_plus_strand: bool

    def gff_range(self) -> str:
        start, end = geenuff_to_gff_start_end(self.start, self.end, self.is_plus_strand)
        return f'{start}-{end}'


@dataclass
class TranscriptRow:
    gene_id: str
    transcript_id: str
    selected: bool
    excluded_from_export: str | None
    seqid: str = ''
    transcript_feature: Feature | None = None
    coding: bool = False
    masks: list[Feature] = field(default_factory=list)
    errors: set[str] = field(default_factory=set)

    @property
    def exported(self) -> bool:
        return self.selected and self.excluded_from_export is None

    @property
    def checked(self) -> bool:
        return self.coding and self.excluded_from_export not in UNCHECKED_REASONS

    @property
    def masked_in_full(self) -> bool:
        return self.transcript_feature is not None and is_masked_whole(self.transcript_feature, self.masks)

    def values(self) -> list[str]:
        tf = self.transcript_feature
        start, end = geenuff_to_gff_start_end(tf.start, tf.end, tf.is_plus_strand) if tf else ('', '')
        return [self.gene_id, self.transcript_id, self.seqid, ('+' if tf.is_plus_strand else '-') if tf else '',
                str(start), str(end), str(self.coding), str(self.selected), str(self.exported),
                self.excluded_from_export or '', ','.join(sorted(self.errors)) if self.checked else 'not checked',
                error_severity(self.errors) if self.checked else '', str(self.masked_in_full),
                ','.join(m.gff_range() for m in self.masks)]


def read_transcripts(con: sqlite3.Connection) -> list[TranscriptRow]:
    """Every transcript of the database with its features, error types and gene."""
    rows = {t_id: TranscriptRow(gene_id, transcript_id, bool(longest), excluded)
            for t_id, transcript_id, longest, gene_id, excluded in con.execute(
                """SELECT t.id, t.given_name, t.longest, sl.given_name, sl.excluded_from_export
                   FROM transcript t JOIN super_locus sl ON sl.id = t.super_locus_id""")}
    for t_id, feature_type, start, end, plus, seqid in con.execute(
            """SELECT tp.transcript_id, f.type, f.start, f."end", f.is_plus_strand, c.seqid
               FROM feature f
               JOIN association_transcript_piece_to_feature a ON a.feature_id = f.id
               JOIN transcript_piece tp ON tp.id = a.transcript_piece_id
               JOIN coordinate c ON c.id = f.coordinate_id
               WHERE f.type IN (?, ?, ?)""", (types.GEENUFF_TRANSCRIPT, types.GEENUFF_CDS, types.GEENUFF_MASK)):
        row, feature = rows[t_id], Feature(start, end, bool(plus))
        if feature_type == types.GEENUFF_TRANSCRIPT:
            row.transcript_feature, row.seqid = feature, seqid
        elif feature_type == types.GEENUFF_CDS:
            row.coding = True
        else:
            row.masks.append(feature)
    for t_id, error_type in con.execute('SELECT transcript_id, type FROM transcript_error'):
        rows[t_id].errors.add(error_type)
    return list(rows.values())


def matches(row: TranscriptRow, args: argparse.Namespace) -> bool:
    if not (row.coding or args.include_noncoding):
        return False
    if args.selected and not row.selected:
        return False
    if args.seqid and row.seqid not in args.seqid:
        return False
    if args.exported is not None and row.exported != args.exported:
        return False
    if args.excluded_reason and row.excluded_from_export not in args.excluded_reason:
        return False
    if args.masked_in_full is not None and row.masked_in_full != args.masked_in_full:
        return False
    if args.error or args.errors is not None or args.severity:
        if not row.checked:
            return False
        if args.error and not row.errors & set(args.error):
            return False
        if args.errors is not None and bool(row.errors) != args.errors:
            return False
        if args.severity and error_severity(row.errors) not in args.severity:
            return False
    return True


def write_counts(rows: list[TranscriptRow], out) -> None:
    """Writes how many of the transcripts and their genes have each error type and severity."""
    by_error, genes_by_error = Counter(), defaultdict(set)
    for row in rows:
        if row.checked:
            for key in [f'error_severity {error_severity(row.errors)}'] + [f'error {e}' for e in row.errors] + [
                    'error-free' if not row.errors else 'with errors']:
                by_error[key] += 1
                genes_by_error[key].add(row.gene_id)
        else:
            by_error['not checked'] += 1
            genes_by_error['not checked'].add(row.gene_id)
    out.write('category\ttranscripts\tgenes\n')
    out.write(f'all listed\t{len(rows)}\t{len({row.gene_id for row in rows})}\n')
    for key in sorted(by_error):
        out.write(f'{key}\t{by_error[key]}\t{len(genes_by_error[key])}\n')


def main(args: argparse.Namespace) -> None:
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')
    logger.info(f'Reading "{args.db_path_in}"')
    con = sqlite3.connect(args.db_path_in)
    rows = [row for row in read_transcripts(con) if matches(row, args)]
    con.close()
    out = open(args.out, 'w') if args.out else sys.stdout
    if args.output == 'table':
        out.write('\t'.join(COLUMNS) + '\n')
        for row in rows:
            out.write('\t'.join(row.values()) + '\n')
    elif args.output in ('gene_ids', 'transcript_ids'):
        ids = [row.gene_id if args.output == 'gene_ids' else row.transcript_id for row in rows]
        out.write(''.join(f'{i}\n' for i in dict.fromkeys(ids)))
    else:
        write_counts(rows, out)
    if args.out:
        out.close()
    logger.info(f'Listed {len(rows)} transcripts of {len({row.gene_id for row in rows})} genes')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--db-path-in', type=str, required=True, help='Path to the GeenuFF SQLite database.')
    parser.add_argument('-o', '--out', type=str, help='output file (default is stdout)')
    parser.add_argument('--output', choices=['table', 'gene_ids', 'transcript_ids', 'counts'], default='table',
                        help=('table: one row per transcript (default); gene_ids or transcript_ids: the IDs alone, '
                              'one per line, e.g. for grep -f; counts: transcripts and genes per error type and '
                              'severity'))
    parser.add_argument('--error', nargs='+', choices=ERROR_TYPES, metavar='TYPE',
                        help=f'transcripts with any of these error types: {", ".join(ERROR_TYPES)}')
    errors = parser.add_mutually_exclusive_group()
    errors.add_argument('--with-errors', dest='errors', action='store_const', const=True,
                        help='transcripts with at least one error')
    errors.add_argument('--error-free', dest='errors', action='store_const', const=False,
                        help='transcripts checked and found without errors')
    parser.add_argument('--severity', nargs='+', choices=SEVERITIES,
                        help=('what the transcript\'s own errors take away: none, flank_masked (only beside it, '
                              'for a missing UTR) or masked_outright'))
    masked = parser.add_mutually_exclusive_group()
    masked.add_argument('--masked-in-full', dest='masked_in_full', action='store_const', const=True,
                        help='transcripts covered entirely by a mask, giving no labels')
    masked.add_argument('--labeled', dest='masked_in_full', action='store_const', const=False,
                        help='transcripts not masked in full, giving labels')
    exported = parser.add_mutually_exclusive_group()
    exported.add_argument('--exported', dest='exported', action='store_const', const=True,
                          help='the transcripts an export takes: one per gene, of genes not excluded')
    exported.add_argument('--not-exported', dest='exported', action='store_const', const=False,
                          help='all other transcripts: other isoforms and those of excluded genes')
    parser.add_argument('--excluded-reason', nargs='+', choices=sorted(UNEXPORTED_REASONS), metavar='REASON',
                        help=('transcripts of genes left out of exports for any of these reasons: '
                              + '; '.join(f'{reason}: {text}' for reason, text in sorted(UNEXPORTED_REASONS.items()))))
    parser.add_argument('--selected', action='store_true',
                        help='only the transcript selected per gene (the longest coding one), not every isoform')
    parser.add_argument('--seqid', nargs='+', help='only transcripts on these sequences')
    parser.add_argument('--include-noncoding', action='store_true',
                        help='list non-coding transcripts too (never checked for errors)')
    args = parser.parse_args()
    require_input_files(parser, [args.db_path_in])
    check_output_paths(parser, [args.out], [args.db_path_in])
    main(args)
