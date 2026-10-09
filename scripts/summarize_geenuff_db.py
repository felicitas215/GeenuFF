#! /usr/bin/env python3
"""Rebuilds as much of the import summary as the database still holds, for when the import log is
gone, in the layout of the import summary. Queries through plain sqlite rather than the ORM, and
fills the counts into geenuff's own ImportStatistics, so the sections and their wording cannot
drift from what the importer logs.

Counts that only existed while the GFF3 was being read are not in the database and cannot be
recovered; the script says which at the end.
"""
import argparse
import logging
import sqlite3
from collections import defaultdict
from types import SimpleNamespace

from geenuff.base import types
from geenuff.applications.importer import (ImportStatistics, GFFErrorHandling, NOT_MASKING_IN_FULL,
                                           exported_outcome, format_summary)

logger = logging.getLogger(__name__)

# an overlap-dropped locus was still part of the sweep at import time, an unplaceable one was not
SWEPT = (types.OVERLAP_DROPPED,)


def has_table(con: sqlite3.Connection, table: str) -> bool:
    return bool(con.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",
                            (table,)).fetchall())


def features_by_transcript(con: sqlite3.Connection, feature_type: str) -> dict[int, list[SimpleNamespace]]:
    """{transcript id: [its features of one type, with start, end, is_plus_strand and coordinate_id]}"""
    out = defaultdict(list)
    for t_id, start, end, plus, coord in con.execute(
            """SELECT tp.transcript_id, f.start, f.end, f.is_plus_strand, f.coordinate_id
               FROM feature f
               JOIN association_transcript_piece_to_feature a ON a.feature_id = f.id
               JOIN transcript_piece tp ON tp.id = a.transcript_piece_id
               WHERE f.type = ?""", (feature_type,)):
        out[t_id].append(SimpleNamespace(start=start, end=end, is_plus_strand=plus, coordinate_id=coord))
    return out


def overlap_partners(spans: dict[int, SimpleNamespace]) -> dict[int, set[int]]:
    """Re-runs the import's overlap sweep over the exported transcripts' spans, {gene id: transcript
    feature}, giving back each gene's partners on the same sequence and strand."""
    by_strand = defaultdict(list)
    for sl_id, tf in spans.items():
        by_strand[(tf.coordinate_id, tf.is_plus_strand)].append((min(tf.start, tf.end), max(tf.start, tf.end),
                                                                sl_id))
    partners = defaultdict(set)
    for strand_spans in by_strand.values():
        for i, j, _, _ in GFFErrorHandling._overlapping_pairs(strand_spans):
            partners[i].add(j)
            partners[j].add(i)
    return partners


def statistics(con: sqlite3.Connection) -> ImportStatistics:
    """The counts of the import summary that the database holds."""
    stats = ImportStatistics()
    stats.total_super_loci = con.execute('SELECT COUNT(*) FROM super_locus').fetchone()[0]
    stats.empty_super_loci = con.execute(
        'SELECT COUNT(*) FROM super_locus sl WHERE NOT EXISTS '
        '(SELECT 1 FROM transcript t WHERE t.super_locus_id = sl.id)').fetchone()[0]
    stats.total_transcripts = con.execute('SELECT COUNT(*) FROM transcript').fetchone()[0]

    logger.info(f'Found {stats.total_super_loci} genes with {stats.total_transcripts} transcripts')

    errors = defaultdict(set)
    for t_id, error_type in con.execute('SELECT transcript_id, type FROM transcript_error'):
        errors[t_id].add(error_type)
    stats.floating_cds_transcripts = sum(types.FLOATING_CDS in e for e in errors.values())

    # {gene id: (exclusion reason, [(coding transcript id, longest)])}
    genes = {}
    for sl_id, reason, t_id, longest in con.execute(
            f"""SELECT sl.id, sl.excluded_from_export, t.id, t.longest FROM transcript t
                JOIN super_locus sl ON sl.id = t.super_locus_id
                WHERE t.id IN (SELECT DISTINCT tp.transcript_id FROM transcript_piece tp
                               JOIN association_transcript_piece_to_feature a ON a.transcript_piece_id = tp.id
                               JOIN feature f ON f.id = a.feature_id WHERE f.type = ?)""",
            (types.GEENUFF_CDS,)):
        genes.setdefault(sl_id, (reason, []))[1].append((t_id, longest))

    logger.info(f'Found {len(genes)} genes with coding transcripts and {len(errors)} transcripts with errors')

    transcript_features = features_by_transcript(con, types.GEENUFF_TRANSCRIPT)
    masks = features_by_transcript(con, types.GEENUFF_MASK)
    selected = {sl_id: next(t_id for t_id, longest in coding if longest)
                for sl_id, (reason, coding) in genes.items() if reason in (None,) + SWEPT}
    logger.info(f'Re-running the overlap sweep over the longest transcripts of {len(selected)} genes')
    partners = overlap_partners({sl_id: transcript_features[t_id][0] for sl_id, t_id in selected.items()})

    for sl_id, (reason, coding) in genes.items():
        if reason is not None:
            stats.unexported_coding_genes[reason] += 1
            stats.unexported_coding_transcripts[reason] += len(coding)
            if reason == types.OVERLAP_DROPPED:
                stats.overlap_loci_dropped_labeled += not errors[selected[sl_id]] - NOT_MASKING_IN_FULL
            continue
        t_id = selected[sl_id]
        stats.longest_transcripts += 1
        stats.unselected_coding_transcripts += len(coding) - 1
        stats.longest_error_free_transcripts += not errors[t_id]
        for error_type in errors[t_id]:
            stats.errors[error_type] += 1
        stats.exported_outcomes[exported_outcome(errors[t_id], transcript_features[t_id][0], masks[t_id],
                                                 len(partners[sl_id]) <= 1)] += 1

    # the pairs of exported genes, settled as at import: in a chain, resolved by dropping one,
    # nested, or both masked outright for their own errors
    for i, j in {tuple(sorted((a, b))) for a, bs in partners.items() for b in bs}:
        spans = [transcript_features[selected[k]][0] for k in (i, j)]
        if len(partners[i]) > 1 or len(partners[j]) > 1:
            stats.overlap_pairs_in_chains += 1
        elif types.OVERLAP_DROPPED in (genes[i][0], genes[j][0]):
            stats.overlap_pairs_resolved += 1
        elif GFFErrorHandling._is_nested(*[sorted((tf.start, tf.end)) for tf in spans]):
            stats.overlap_pairs_nested += 1
        else:
            stats.overlap_pairs_both_masked_outright += 1

    stats.overlap_pairs_recorded = con.execute('SELECT COUNT(*) FROM super_locus_overlap').fetchone()[0]
    stats.super_loci_in_overlap_pairs = con.execute(
        """SELECT COUNT(*) FROM (SELECT super_locus_id FROM super_locus_overlap
                                 UNION SELECT partner_id FROM super_locus_overlap)""").fetchone()[0]
    return stats


def main(args: argparse.Namespace) -> None:
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')
    logger.info(f'Reading "{args.db_path_in}"')
    con = sqlite3.connect(args.db_path_in)
    if not has_table(con, 'transcript_error'):
        raise SystemExit(f'"{args.db_path_in}" was written by a GeenuFF without the transcript_error '
                         f'table; import it again to summarize it')
    species = con.execute('SELECT species FROM genome').fetchone()
    stats = statistics(con)
    con.close()

    title = f'summary rebuilt from "{args.db_path_in}"' + (f' for "{species[0]}"' if species else '')
    print(format_summary(title, stats.summary_sections(reading_counts=False), stats.summary_notes()), flush=True)
    print('  not in the database, only ever in the import log:', flush=True)
    for line in ('GFF lines skipped or left out while grouping them into genes, and the coding '
                 'transcripts lost through them',
                 'genes created for transcripts naming a Parent ID that matches no line or to mask '
                 'floating CDS lines, and gene lines on neither strand, indistinguishable from other '
                 'genes',
                 'transcripts with exons built from their CDS and UTR lines, or with a stop codon added',
                 'transcripts dropped, never written at all'):
        print(f'    - {line}', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--db-path-in', type=str, required=True,
                        help='Path to the GeenuFF SQLite database.')
    main(parser.parse_args())
