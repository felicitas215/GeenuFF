#! /usr/bin/env python3
"""Rebuilds as much of the import summary as the database still holds, for when the import log is
gone. Queries through plain sqlite rather than the ORM, so a database an older GeenuFF wrote still
opens and the sections its schema cannot answer are simply left out. The feature types and
exclusion reasons come from geenuff itself, so they cannot drift from what the importer writes.

Counts that only existed while the GFF3 was being read are not in the database and cannot be
recovered; the script says which at the end.
"""
import argparse
import sqlite3
from collections import defaultdict

from geenuff.base import types
from geenuff.applications.importer import UNEXPORTED_REASONS

ERROR_TYPES = tuple(types.geenuff_error_type_values)
TRANSCRIPT_TYPE = types.GEENUFF_TRANSCRIPT
CDS_TYPE = types.GEENUFF_CDS
# an overlap-dropped locus was still part of the sweep at import time, an unplaceable one was not
SWEPT = (types.OVERLAP_DROPPED,)

# a locus id and exclusion reason, its transcript count, and the span of its exported transcript
LocusRow = tuple[int, str | None, int, int | None, int | None, bool | None, int | None]
# a counted line of the summary, or a heading standing over the lines that follow it
Entry = tuple[int, str] | str


def has_column(con: sqlite3.Connection, table: str, column: str) -> bool:
    return any(row[1] == column for row in con.execute(f'PRAGMA table_info({table})'))


def has_table(con: sqlite3.Connection, table: str) -> bool:
    return bool(con.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
                            (table,)).fetchall())


def locus_rows(con: sqlite3.Connection, excluded_column: bool) -> list[LocusRow]:
    """Per super locus: its exclusion reason, transcript count, and the span, coordinate and
    strand of the one transcript an export emits for it, where it has one."""
    reason_col = 'sl.excluded_from_export' if excluded_column else 'NULL'
    query = f"""SELECT sl.id, {reason_col},
                       (SELECT COUNT(*) FROM transcript t WHERE t.super_locus_id = sl.id),
                       f.start, f.end, f.is_plus_strand, f.coordinate_id
                FROM super_locus sl
                LEFT JOIN transcript t ON t.super_locus_id = sl.id AND t.longest = 1
                LEFT JOIN transcript_piece tp ON tp.transcript_id = t.id
                LEFT JOIN association_transcript_piece_to_feature a
                    ON a.transcript_piece_id = tp.id
                LEFT JOIN feature f ON f.id = a.feature_id AND f.type = ?
                GROUP BY sl.id"""
    return list(con.execute(query, (TRANSCRIPT_TYPE,)))


def overlap_sweep(loci: list[LocusRow]) -> dict[int, set[int]]:
    """Re-runs the import's overlap sweep over the spans the database holds, giving back each
    locus' partners. Loci an export never reached are left out, as they were at import."""
    groups = defaultdict(list)
    for sl_id, reason, _, start, end, plus, coord in loci:
        if start is None or (reason is not None and reason not in SWEPT):
            continue
        groups[(coord, plus)].append((min(start, end), max(start, end), sl_id))
    partners = defaultdict(set)
    for spans in groups.values():
        active = []
        for lo, hi, sl_id in sorted(spans):
            active = [(a_hi, a_id) for a_hi, a_id in active if a_hi > lo]
            for _, a_id in active:
                partners[sl_id].add(a_id)
                partners[a_id].add(sl_id)
            active.append((hi, sl_id))
    return partners


def report(con: sqlite3.Connection, out: list[Entry]) -> None:
    """Appends every section the database can answer to out."""
    excluded_column = has_column(con, 'super_locus', 'excluded_from_export')
    loci = locus_rows(con, excluded_column)
    excluded = defaultdict(int)
    n_empty = 0
    for _, reason, n_transcripts, *_ in loci:
        if reason is not None:
            excluded[reason] += 1
        if not n_transcripts:
            n_empty += 1
    distinct_names = con.execute('SELECT COUNT(DISTINCT given_name) FROM super_locus').fetchone()[0]

    out.append('loci:')
    out.append((len(loci), 'gene lines in the database'))
    out.append((n_empty, 'of them with no transcript under them'))
    out.append((len(loci) - distinct_names, 'gene IDs used by more than one gene line'))
    for reason, count in sorted(excluded.items()):
        out.append((count, f'genes kept but never exported: '
                           f'{UNEXPORTED_REASONS.get(reason, reason)}'))

    coding = {r[0] for r in con.execute(
        """SELECT DISTINCT t.id FROM transcript t
           JOIN transcript_piece tp ON tp.transcript_id = t.id
           JOIN association_transcript_piece_to_feature a ON a.transcript_piece_id = tp.id
           JOIN feature f ON f.id = a.feature_id WHERE f.type = ?""", (CDS_TYPE,))}
    n_transcripts = con.execute('SELECT COUNT(*) FROM transcript').fetchone()[0]
    placeholders = ','.join('?' * len(ERROR_TYPES))
    exported = list(con.execute(
        f"""SELECT t.id, COUNT(f.id)
            FROM transcript t
            JOIN super_locus sl ON sl.id = t.super_locus_id
            LEFT JOIN transcript_piece tp ON tp.transcript_id = t.id
            LEFT JOIN association_transcript_piece_to_feature a
                ON a.transcript_piece_id = tp.id
            LEFT JOIN feature f ON f.id = a.feature_id AND f.type IN ({placeholders})
            WHERE t.longest = 1
              {'AND sl.excluded_from_export IS NULL' if excluded_column else ''}
            GROUP BY t.id""", ERROR_TYPES))

    out.append('transcripts:')
    out.append((n_transcripts, 'transcripts in the database'))
    out.append((len(coding), 'of them with a CDS'))
    out.append((len(exported), 'selected for export, one per locus that will be exported'))
    out.append((sum(1 for _, n_errors in exported if not n_errors),
                'of those selected with no error at all'))

    if has_table(con, 'super_locus_overlap'):
        pairs = con.execute('SELECT COUNT(*) FROM super_locus_overlap').fetchone()[0]
        in_pairs = con.execute("""SELECT COUNT(*) FROM (
                                    SELECT super_locus_id AS id FROM super_locus_overlap
                                    UNION SELECT partner_id FROM super_locus_overlap)""").fetchone()[0]
        out.append('overlapping coding genes, over every coding transcript they have:')
        out.append((in_pairs, f'genes recorded in the super_locus_overlap table, in {pairs} pairs'))

    partners = overlap_sweep(loci)
    exported_pairs = {tuple(sorted((a, b))) for a, bs in partners.items() for b in bs}
    chained = {p for p in exported_pairs if len(partners[p[0]]) > 1 or len(partners[p[1]]) > 1}
    n_resolved = excluded.get(types.OVERLAP_DROPPED, 0)
    out.append('overlapping genes, over the one transcript each of them exports:')
    out.append((len(partners), f'genes overlap at least one other, in {len(exported_pairs)} pairs'))
    out.append((n_resolved, 'pairs where one gene was kept whole and its partner dropped'))
    out.append((len(exported_pairs) - len(chained) - n_resolved,
                'pairs where neither gene could be kept, so both are masked'))
    out.append((len(chained), 'pairs left alone, one of the two overlapping a further gene'))

    errors = list(con.execute(
        f"""SELECT f.type, COUNT(DISTINCT t.id) FROM feature f
            JOIN association_transcript_piece_to_feature a ON a.feature_id = f.id
            JOIN transcript_piece tp ON tp.id = a.transcript_piece_id
            JOIN transcript t ON t.id = tp.transcript_id
            WHERE f.type IN ({placeholders}) GROUP BY f.type ORDER BY f.type""", ERROR_TYPES))
    if errors:
        out.append('errors, counted once per transcript they occur in:')
        out += [(count, error_type) for error_type, count in errors]


def main(args: argparse.Namespace) -> None:
    con = sqlite3.connect(args.db_path_in)
    species = con.execute('SELECT species FROM genome').fetchone()
    out: list[Entry] = [f'summary rebuilt from "{args.db_path_in}"'
                        + (f' for "{species[0]}"' if species else '') + ':']
    report(con, out)
    con.close()

    width = max((len(str(entry[0])) for entry in out if isinstance(entry, tuple)), default=1)
    for entry in out:
        if isinstance(entry, tuple):
            print(f'    {entry[0]:>{width}}  {entry[1]}')
        else:
            print(entry if entry.startswith('summary') else f'  {entry}')

    print('\n  not in the database, only ever in the import log:')
    for line in ('gene lines on neither strand, saved without transcripts and so '
                 'indistinguishable from any other empty gene',
                 'transcripts dropped for having no storable range, never written at all',
                 'exon/CDS lines dropped while reading the GFF3, whether as duplicates, as '
                 'belonging to a transcript-less gene, or as an isoform with no transcript line',
                 'lines skipped for a feature type GeenuFF has no use for'):
        print(f'    - {line}')
    print('\n  errors are counted from the features stored, so one whose masked range came out\n'
          '  empty still counts, but one discarded for running backwards does not')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--db-path-in', type=str, required=True,
                        help='Path to the GeenuFF SQLite database.')
    main(parser.parse_args())
