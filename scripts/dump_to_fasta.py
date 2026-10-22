#! /usr/bin/env python3
import logging

from geenuff.applications.exporters.sequence import FastaExportController
from geenuff.applications.exporter import RangeArgParser, MODES


def main(args):
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')
    controller = FastaExportController(args.db_path_in, args.longest, args.coding)
    if args.mode in MODES:
        controller.prep_ranges(MODES[args.mode])
    else:
        raise NotImplementedError(f'Requested mode ({args.mode}) not in implemented types {list(MODES.keys())}')

    controller.write_fa(args.out)


if __name__ == '__main__':
    parser = RangeArgParser()
    parser.parse_args()
    main(parser.args)
