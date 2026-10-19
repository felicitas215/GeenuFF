#! /usr/bin/env python3
import argparse

from geenuff.applications.exporters.gff3 import FilteredGff3ExportController
from geenuff.base.helpers import check_output_paths, require_input_files, setup_logging


def main(args: argparse.Namespace) -> None:
    setup_logging(args.log_file)
    controller = FilteredGff3ExportController(args.db_path_in)
    controller.write_filtered_gff3(args.out, include_erroneous=args.include_erroneous)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=('Writes a GFF3 file containing one transcript per gene, the longest coding '
                     'one, for comparing a Helixer prediction against its reference. By default '
                     'only the ones not masked whole, i.e. the same gene models Helixer\'s h5 '
                     'export gives labels from.'))
    parser.add_argument('--db-path-in', type=str, required=True,
                        help='Path to the GeenuFF SQLite input database.')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='output GFF3 file path')
    parser.add_argument('--include-erroneous', action='store_true',
                        help=('write every gene that can be written at all, whatever is wrong '
                              'with it, including ones dropped from the h5 export for '
                              'overlapping another gene. Only genes whose features cannot be '
                              'placed are left out, and the log says how many'))
    parser.add_argument('--log-file', type=str, help='file to write the log to, in addition to the terminal')
    args = parser.parse_args()
    require_input_files(parser, [args.db_path_in])
    check_output_paths(parser, [args.out, args.log_file], [args.db_path_in])
    main(args)
