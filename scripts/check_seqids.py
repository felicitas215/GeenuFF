#! /usr/bin/env python3
"""Checks whether the sequence IDs of a GFF3 file are those of its FASTA file (the header up to the
first space). Where they are not, GeenuFF's import and many other tools fail or leave the genes on
those sequences out; the log says what to check. Exits with 1 where IDs do not match, 0 otherwise.
"""
import argparse
import logging
import os

from geenuff.base.helpers import fasta_ids, get_seqids_from_gff

logger = logging.getLogger(__name__)

MAX_LISTED = 5


def listed(ids: set[str]) -> str:
    shown = sorted(ids)[:MAX_LISTED]
    return ', '.join(shown) + (', ...' if len(ids) > MAX_LISTED else '')


def main(args: argparse.Namespace) -> None:
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')
    logger.info(f'Checking the sequence IDs of "{os.path.basename(args.gff)}" against '
                f'"{os.path.basename(args.fasta)}"')
    fasta, gff = set(fasta_ids(args.fasta)), get_seqids_from_gff(args.gff)
    missing = gff - fasta
    if not missing:
        logger.info(f'All {len(gff)} sequence IDs of the GFF3 file are FASTA IDs: the files match by ID '
                    f'({len(fasta - gff)} FASTA sequences have no GFF3 lines)')
        return
    logger.warning(f'{len(missing)} of {len(gff)} sequence IDs of the GFF3 file are no FASTA ID, e.g. '
                   f'{listed(missing)} (FASTA IDs e.g. {listed(fasta)}). GeenuFF\'s import and many other tools '
                   f'fail or leave out the genes on them. Check that the FASTA and GFF3 files come from the '
                   f'same source and assembly; if the sequences were renamed, e.g. on upload to NCBI, '
                   f'scripts/rename_gff_seqids.py pairs and renames them.')
    raise SystemExit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--fasta', type=str, required=True, help='FASTA file (may be gzipped)')
    parser.add_argument('--gff', type=str, required=True, help='GFF3 file (may be gzipped)')
    main(parser.parse_args())
