#!/usr/bin/env python
import sys
import argparse
import logging
_logger = logging.getLogger(__name__)
import random


def parse_args(args):
    """Parse command line parameters"""
    parser = argparse.ArgumentParser(
        description="Parse parameters for %s" % __file__)
    parser.add_argument(
        "genomesize",
        type=int,
        help="Genome or chromosome length",
        metavar="genome_size")
    parser.add_argument(
        "genomename",
        type=str,
        help="Genome or chromosome name",
        metavar="genome_name")
    parser.add_argument(
        "outfile",
        type=str,
        help="Output junction detailed csv file from rnatk_ja",
        metavar="out_jd.csv")
    parser.add_argument(
        "-n",
        "--num",
        type=int,
        dest="num",
        default=1000)
    parser.add_argument(
        "-g",
        "--mingap",
        type=int,
        dest="mingap",
        help="Minimum length of junction gap",
        default=30)
    parser.add_argument(
        "-a",
        "--armsize",
        type=int,
        dest="armsize",
        help="Minimum size of flanking arms",
        default=20)
    parser.add_argument(
        "-s",
        "--seed",
        type=int,
        dest="seed",
        help="Seed for random number generation",
        default=2020)
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
    random.seed(args.seed)
    a = 0 + args.armsize
    b = args.genomesize - args.armsize
    b = args.genomesize - args.armsize
    with open(args.outfile, "w") as foh:
        k = 0
        while k < args.num:
            i = random.randint(a, b)
            j = random.randint(a, b)
            if i > j:
                i, j = j, i

            if j - i < args.mingap:
                continue

            k += 1
            foh.write("\t".join(map(str, [
                k, 1, args.genomename, i, j, "%d|1" % (k, ), 1, '+',
                i, j, '0,0,0', 2, "1,1,", "0,%d," % (j-i-1, ),
                1, 1
            ]))+"\n")

    _logger.info("Ends.")


def run():
    """Entry point for console_scripts"""
    main(sys.argv[1:])


if __name__ == '__main__':
    run()
