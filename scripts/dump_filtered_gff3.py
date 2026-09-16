#! /usr/bin/env python3
import argparse
import logging

from geenuff.applications.exporters.gff3 import FilteredGff3ExportController


def main(args: argparse.Namespace) -> None:
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format='%(asctime)s - %(levelname)s: %(message)s',
                        datefmt='%d-%b-%y %H:%M:%S')
    controller = FilteredGff3ExportController(args.db_path_in)
    controller.write_filtered_gff3(args.out)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=('Writes a GFF3 file containing only the longest, fully error-free '
                     'transcript per super locus. These are the same gene models Helixer\'s '
                     'h5 export uses, i.e. the ones used for training.'))
    parser.add_argument('--db-path-in', type=str, required=True,
                        help='Path to the GeenuFF SQLite input database.')
    parser.add_argument('-o', '--out', type=str,
                        help='output GFF3 file path (default is stdout)')
    parser.add_argument('--verbose', action='store_true',
                        help='log at DEBUG level instead of INFO')

    main(parser.parse_args())
