#!/usr/bin/env python3
import os
import sys
import yaml
import logging
import argparse

from geenuff.applications.importer import ImportController
from geenuff.base.helpers import check_output_paths, require_input_files

logger = logging.getLogger(__name__)

DEFAULT_CONFIG_FILE = 'config/import.yml'


class PathFinder(object):
    INPUT = 'input'
    OUTPUT = 'output'

    def __init__(self, db_path, basedir, species, fasta=None, gff=None, logfile=None):
        # directories
        self.db_out = db_path
        self.basedir = basedir
        self.species = species
        self.input = '{}/{}/'.format(self.basedir, PathFinder.INPUT)
        self.output = '{}/{}/'.format(self.basedir, PathFinder.OUTPUT)
        # basedir and, unless both input files are given, basedir/input exist (checked before)
        if self.basedir is not None:
            os.makedirs(self.output, exist_ok=True)
        # files
        self.fasta_in = self._get_fa(fasta)
        self.gff_in = self._get_gff(gff)
        if not self.db_out:
            self.db_out = '{}{}.sqlite3'.format(self.output, species)
        # beside the database, named after it, e.g. genome.sqlite3 -> genome_geenuff_import.log
        db_stem = self.db_out.removeprefix('sqlite:///').removesuffix('.sqlite3')
        self.problems_out = logfile or f'{db_stem}_geenuff_import.log'

    def _get_fa(self, provided):
        if provided is not None:
            return provided
        maybe = os.listdir(self.input)
        # todo, actual file type detection
        maybe = [x for x in maybe if any([x.endswith(ok) for ok in ['.fa', '.fasta', '.fa.gz', '.fasta.gz',
                                                                    '.fa.zip', '.fasta.zip']])]
        self._confirm_exactly_one(maybe, 'fasta')
        return self.input + maybe[0]

    def _get_gff(self, provided):
        if provided is not None:
            return provided
        maybe = os.listdir(self.input)
        maybe = [x for x in maybe if any([x.endswith(ok) for ok in ['.gff', '.gff3', '.gff.gz', '.gff3.gz',
                                                                    '.gff.zip', '.gff3.zip']])]
        self._confirm_exactly_one(maybe, 'gff')
        return self.input + maybe[0]

    @staticmethod
    def _confirm_exactly_one(possibilities, info):
        assert len(possibilities) == 1, 'no(n) unique {} file found as input. Found: {}'.format(
                                            info, possibilities)


def main(args, paths):
    # checked before the log file is opened, which would overwrite the existing database's import log
    if os.path.exists(paths.db_out.removeprefix('sqlite:///')) and not args.replace_db:
        sys.exit(f'database already existing at {paths.db_out} and --replace-db not set')

    msg_fmt_str = '%(asctime)s - %(levelname)s: %(message)s'
    logging.basicConfig(filename=paths.problems_out,
                        filemode='w',
                        level=logging.DEBUG if args.verbose else logging.INFO,
                        format=msg_fmt_str)
    # log to file and stderr simultaneously
    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setFormatter(logging.Formatter(fmt=msg_fmt_str))
    logging.getLogger().addHandler(stdout_handler)

    # try to load config
    config = {}
    if os.path.isfile(args.config_file):
        with open(args.config_file, 'r') as f:
            try:
                # an empty file loads as None
                config = yaml.safe_load(f) or {}
            except yaml.YAMLError as e:
                logger.error(f'An error occurred during parsing of the YAML config file: {e}')
                sys.exit(1)
    elif args.config_file != DEFAULT_CONFIG_FILE:
        logger.error(f'Config file "{args.config_file}" not found')
        sys.exit(1)
    else:
        logger.info('No config file found, using default values')

    controller = ImportController(database_path=paths.db_out, config=config, replace_db=args.replace_db)
    genome_args = {}
    for key in ['species', 'accession', 'version', 'acquired_from']:
        genome_args[key] = vars(args)[key]
    controller.add_genome(paths.fasta_in, paths.gff_in, genome_args)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--basedir', help='organized output (& input) directory. If this is not set, --fasta, --gff3 '
                        'and --db-path must be set.')
    parser.add_argument('--config-file', type=str, default=DEFAULT_CONFIG_FILE,
                        help=f'YAML file setting min_intron_length and min_cds_length, the intron and spliced CDS '
                             f'lengths in bp below which a transcript is masked, e.g. config/import.yml of the '
                             f'GeenuFF repository (default: {DEFAULT_CONFIG_FILE} in the working directory if '
                             f'present, else 20 and 60)')

    custominput = parser.add_argument_group('Override default with custom input/output location')
    custominput.add_argument('--gff3', help='gff3 formatted file to parse / standardize')
    custominput.add_argument('--fasta', help='fasta file to parse / standardize')
    custominput.add_argument('--db-path', help=('path of the GeenuFF database, which must end in .sqlite3 (default: '
                                                'basedir/output/<species>.sqlite3)'))
    custominput.add_argument('--log-file', help=('output path for import log (default: beside the database, named '
                                                 'after it, e.g. genome.sqlite3 -> genome_geenuff_import.log)'))

    parser.add_argument('--replace-db', action='store_true',
                        help=('whether to override a GeenuFF database found at '
                              'the default location or at the location of --db-path'))
    parser.add_argument('--verbose', action='store_true',
                        help=('log every individual error as it is found (DEBUG level), '
                              'instead of just the aggregate summary at the end, and the progress lines of '
                              'every sequence, of which only 30 are logged otherwise'))

    genome_attr = parser.add_argument_group('Possible genome attributes:')
    genome_attr.add_argument('--species', required=True, help='name of the species')
    genome_attr.add_argument('--accession', default='', help='')
    genome_attr.add_argument('--version', default='', help='genome version')
    genome_attr.add_argument('--acquired-from', default='', help='genome source')

    args = parser.parse_args()
    if args.basedir is None and None in (args.fasta, args.gff3, args.db_path):
        parser.error('without --basedir, --fasta, --gff3 and --db-path are needed')
    if args.basedir is not None:
        # the input files are looked for in basedir/input, unless both are given
        needed = args.basedir if args.fasta and args.gff3 else os.path.join(args.basedir, PathFinder.INPUT)
        if not os.path.isdir(needed):
            parser.error(f'the directory "{needed}" does not exist')
    if args.db_path and not args.db_path.endswith('.sqlite3'):
        parser.error(f'--db-path must end in .sqlite3, got "{args.db_path}"')
    paths = PathFinder(args.db_path, args.basedir, args.species, fasta=args.fasta, gff=args.gff3,
                       logfile=args.log_file)
    require_input_files(parser, [paths.fasta_in, paths.gff_in])
    check_output_paths(parser, [paths.db_out, paths.problems_out], [paths.fasta_in, paths.gff_in])
    main(args, paths)
