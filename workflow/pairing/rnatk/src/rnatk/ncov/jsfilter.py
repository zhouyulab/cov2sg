#!/usr/bin/env python
import sys
import argparse
import logging
_logger = logging.getLogger(__name__)


def parse_args(args):
    """Parse command line parameters"""
    parser = argparse.ArgumentParser(
        description="Parse parameters for %s" % __file__)
    parser.add_argument(
        "jsbed",
        type=str,
        help="BED file for junctions sites",
        metavar="js.bed")
    parser.add_argument(
        "outfile",
        type=str,
        help="Enabled junctions view for IGV",
        metavar="js_igv.bed")
    parser.add_argument(
        "-t",
        "--top",
        type=int,
        dest="top",
        help="Top N junctions",
        default=20)
    parser.add_argument(
        "-g",
        "--mingap",
        type=int,
        dest="mingap",
        help="Minimum length of junction gap",
        default=30)
    parser.add_argument(
        "-s",
        "--strand",
        type=str,
        dest="strand",
        choices=['.', '+', '-'],
        help="strand .|+|-",
        default=".")
    parser.add_argument(
        "-n",
        "--numevi",
        type=int,
        dest="numevi")
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
    _logger.debug("Starting reading ...")
    with open(args.outfile, "w") as foh:
        foh.write("track graphType=junctions\n")
        i = 0
        for line in open(args.jsbed):
            b = line.rstrip().split("\t")
            if args.numevi is not None:
                b[4] = args.numevi
            else:
                b[4] = int(b[12])  # numevi column
            if int(b[2]) - int(b[1]) < args.mingap:
                continue

            if args.strand != '.':
                if b[5] != args.strand:
                    continue

            i += 1
            foh.write("\t".join(map(str, b[:12]))+"\n")
            if i >= args.top:
                break

    _logger.info("Ends.")


def run():
    """Entry point for console_scripts"""
    main(sys.argv[1:])


if __name__ == '__main__':
    run()


