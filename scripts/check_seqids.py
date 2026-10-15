#! /usr/bin/env python3
"""Checks whether the sequence IDs of a GFF3 file are those of its FASTA file (the header up to the
first whitespace). A GFF3 ID contained in exactly one FASTA ID (e.g. Chr1 in Chr1.v2) is paired with
it by GeenuFF's import, with a warning to look at the pairs. Otherwise, GeenuFF's import and many
other tools fail or leave the genes on those sequences out; the log says what to check. Exits with 1
where IDs cannot be paired, 0 otherwise.
"""
import argparse
import logging
import os

from geenuff.base.helpers import NonMatchableIDs, fasta_ids, get_seqids_from_gff, match_seqids, setup_logging

logger = logging.getLogger(__name__)

MAX_LISTED = 5


def listed(ids: set[str]) -> str:
    shown = sorted(ids)[:MAX_LISTED]
    return ', '.join(shown) + (', ...' if len(ids) > MAX_LISTED else '')


def main(args: argparse.Namespace) -> None:
    setup_logging()
    logger.info(f'Checking the sequence IDs of "{os.path.basename(args.gff)}" against '
                f'"{os.path.basename(args.fasta)}"')
    fasta, gff = set(fasta_ids(args.fasta)), get_seqids_from_gff(args.gff)
    missing = gff - fasta
    if not missing:
        logger.info(f'All {len(gff)} sequence IDs of the GFF3 file are FASTA IDs: the files match by ID '
                    f'({len(fasta - gff)} FASTA sequences have no GFF3 lines)')
        return
    try:
        # pairs as GeenuFF's import does, warning with the pairs where they are no exact match
        match_seqids(fasta, gff)
    except NonMatchableIDs:
        logger.warning(f'{len(missing)} of {len(gff)} sequence IDs of the GFF3 file are no FASTA ID, e.g. '
                       f'{listed(missing)} (FASTA IDs e.g. {listed(fasta)}). GeenuFF\'s import and many other '
                       f'tools fail or leave out the genes on them. Check that the FASTA and GFF3 files come from '
                       f'the same source and assembly; if the sequences were renamed, e.g. on upload to NCBI, '
                       f'rename_gff_seqids.py pairs and renames them.')
        raise SystemExit(1)
    logger.info(f'GeenuFF\'s import pairs these {len(missing)} GFF3 sequence IDs with the FASTA IDs containing '
                f'them, but give the pairs a look since they are no exact match')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--fasta', type=str, required=True, help='FASTA file (may be gzipped)')
    parser.add_argument('--gff', type=str, required=True, help='GFF3 file (may be gzipped)')
    main(parser.parse_args())
