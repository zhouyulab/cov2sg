#!/usr/bin/env python
"""
Filter CSV file on one column
"""
import sys
import argparse
from bx.tabular.io import TableReader, Comment
import logging

_logger = logging.getLogger(__name__)
from rnatk import __version__


def csv2filter(fcsv, relation, threshold, col=5, outfile=None):
    """Filter CSV file on given col with relation threshold comparison"""
    is_float = True
    try:
        threshold = float(threshold)
    except ValueError:
        is_float = False

    col = col - 1
    if not outfile:
        foh = sys.stdout
    else:
        foh = open(outfile, "w")

    for row in TableReader(open(fcsv),
                           return_header=False, return_comments=True, comment_lines_startswith=["#"]):
        if isinstance(row, Comment):
            foh.write(str(row) + '\n')
            continue

        try:
            score = row.fields[col]
        except IndexError:
            sys.stderr.write("col index (1-based) is out of range\n")
            return False

        if is_float:
            try:
                score = float(score)
            except ValueError:
                sys.stderr.write('threshold is numeric but not col\n')
                return False
        try:
            if eval('score %s threshold' % relation):
                foh.write(str(row) + '\n')
        except Exception as err:
            sys.stderr.write("Handling run-time error:", err)
            return False

    if not outfile:
        foh.close()
    return True


def parse_args(args):
    """Parse command line parameters"""
    parser = argparse.ArgumentParser(
        description="Parse parameters for rnatk.bin.csv2filter")
    parser.add_argument(
        "fcsv",
        type=str,
        help="Input tab-split CSV file",
        metavar="in.csv")
    parser.add_argument(
        "relation",
        type=str,
        choices=['<', '>', '==', '>=', '<=', '!='],
        help="Relation choice to compare with given column")
    parser.add_argument(
        "threshold",
        type=str,
        help="Threshold for the filtering (float or string)")
    parser.add_argument(
        "-c",
        "--column",
        dest="col",
        type=int,
        default=5,
        help="Column (1-based) to filter (default 5)")
    parser.add_argument(
        "-o",
        "--outfile",
        dest="outfile",
        type=str,
        help="Output file (default: stdout)")
    parser.add_argument(
        '--version',
        action='version',
        version='rnatk {ver}'.format(ver=__version__))
    parser.add_argument(
        '-v',
        '--verbose',
        dest="loglevel",
        help="set loglevel to INFO",
        action='store_const',
        const=logging.INFO)
    parser.add_argument(
        '-vv',
        '--very-verbose',
        dest="loglevel",
        help="set loglevel to DEBUG",
        action='store_const',
        const=logging.DEBUG)
    return parser.parse_args(args)


def setup_logging(loglevel):
    """Setup basic logging
    Args:
      loglevel (int): minimum loglevel for emitting messages
    """
    logformat = "[%(asctime)s] %(levelname)s:%(name)s:%(message)s"
    logging.basicConfig(level=loglevel, stream=sys.stdout,
                        format=logformat, datefmt="%Y-%m-%d %H:%M:%S")


def main(args):
    """Main entry point allowing external calls"""
    args = parse_args(args)
    setup_logging(args.loglevel)
    _logger.debug("Starting reading CSV ...")

    csv2filter(args.fcsv, args.relation, args.threshold, col=args.col, outfile=args.outfile)

    _logger.info("Ends.")


def run():
    """Entry point for console_scripts"""
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
