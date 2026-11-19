#! /usr/bin/env python3
"""Checks whether the sequence IDs of a GFF3 file are those of its FASTA file (the header up to the
first whitespace), as GeenuFF's import requires. Otherwise, GeenuFF's import and many other tools
fail or leave the genes on those sequences out; the log says what to check. Exits with 1 where a
GFF3 ID is no FASTA ID, 0 otherwise.
"""
import argparse
import logging
import os

from geenuff.base.helpers import (add_version_arg, listed, require_fasta_and_gff, require_input_files, run_main,
                                  setup_logging)
from geenuff.parsing.fastahelper import FastaParser
from geenuff.parsing.gffhelper import GFFParser

logger = logging.getLogger(__name__)


def main(args: argparse.Namespace) -> None:
    setup_logging()
    logger.info(f'Checking the sequence IDs of "{os.path.basename(args.gff)}" against '
                f'"{os.path.basename(args.fasta)}"')
    fasta, gff = set(FastaParser.ids(args.fasta)), GFFParser.seqids(args.gff)
    missing = gff - fasta
    if not missing:
        logger.info(f'All {len(gff)} sequence IDs of the GFF3 file are FASTA IDs: the files match by ID '
                    f'({len(fasta - gff)} FASTA sequences have no GFF3 lines)')
        return
    logger.warning(f'{len(missing)} of {len(gff)} sequence IDs of the GFF3 file are no FASTA ID, e.g. '
                   f'{listed(sorted(missing))} (FASTA IDs e.g. {listed(sorted(fasta))}). GeenuFF\'s import '
                   f'and many other tools fail or leave out the genes on them. Check that the FASTA and GFF3 '
                   f'files come from the same source and assembly; if the sequences were renamed, e.g. on '
                   f'upload to NCBI, rename_gff_seqids.py pairs and renames them.')
    raise SystemExit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_version_arg(parser)
    parser.add_argument('--fasta', type=str, required=True, help='FASTA file (may be gzipped)')
    parser.add_argument('--gff', type=str, required=True, help='GFF3 file (may be gzipped)')
    args = parser.parse_args()
    require_input_files(parser, [args.fasta, args.gff])
    require_fasta_and_gff(parser, args.fasta, args.gff)
    run_main(main, args)
