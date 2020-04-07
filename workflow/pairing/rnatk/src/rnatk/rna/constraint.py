#!/usr/bin/env python
"""Add structure constraint in RNAcofold"""
import os
from Bio import SeqIO
import sys
import argparse
import logging
_logger = logging.getLogger(__name__)


BPs = []
for b1, b2 in [('A', 'T'), ('G', 'C'), ('G', 'T'), ('A', 'U'), ('G', 'U')]:
    BPs.extend([(b1, b2), (b2, b1)])
BPs = set(BPs)


def is_paired(b1, b2):
    if (b1, b2) in BPs:
        return True
    return False


def add_constraint(dna, prime=5, struct=None):
    assert prime in (5, 3)
    dna = dna.upper()
    s1, s2 = dna.split("&")
    c1 = ["."] * len(s1)
    c2 = ["."] * len(s2)
    if struct is not None:
        c1, c2 = struct.split("&")
        c1 = list(c1)
        c2 = list(c2)
    for k in range(min(len(s1), len(s2))):
        if prime == 5:
            i, j = k, -(k+1)
        else:
            i, j = -(k+1), k
        matched = False
        if c1[i] == "." and c2[j] == ".":
            if is_paired(s1[i], s2[j]):
                c1[i], c2[j] = '(', ')'
                matched = True
            else:
                break
        if not matched:
            break

    return "".join(c1 + ['&'] + c2)


def add_interconst(dna):
    s1, s2 = dna.split("&")
    return "<"*len(s1) + "&" + ">"*len(s2)


def parse_args(args):
    """Parse command line parameters"""
    parser = argparse.ArgumentParser(
        description="Parse parameters for %s" % __file__)
    parser.add_argument(
        "fa_in",
        type=str,
        help="Input Fasta file",
        metavar="input.fa")
    parser.add_argument(
        "fa_out",
        type=str,
        help="Output Fasta file",
        metavar="out.fa")
    parser.add_argument(
        "-d",
        "--direction",
        type=str,
        dest="direction",
        choices=['lft', 'rgt'],
        help="Side of sequence relative to junc start/end",
        default="rgt")
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

    with open(args.fa_out, "w") as foh:
        for r in SeqIO.parse(args.fa_in, "fasta"):
            s = str(r.seq)
            if args.direction == 'rgt':
                c = add_constraint(s, prime=5)
            else:
                c = add_constraint(s, prime=3)
            foh.write(">%s\n%s\n%s\n" % (r.id, s, c))

    _logger.info("Ends.")


def run():
    """Entry point for console_scripts"""
    main(sys.argv[1:])


def parse_args_cofold(args):
    """Parse command line parameters"""
    parser = argparse.ArgumentParser(
        description="Parse parameters for %s" % __file__)
    parser.add_argument(
        "cofold_in",
        type=str,
        help="Input file in RNAcofold format",
        metavar="cofold.csv")
    parser.add_argument(
        "cofold_out",
        type=str,
        help="Output file in RNAcofold format",
        metavar="out.csv")
    parser.add_argument(
        "-d",
        "--direction",
        type=str,
        dest="direction",
        choices=['lft', 'rgt'],
        help="Side of sequence relative to junc start/end",
        default="rgt")
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


def fix_cofold():
    args = parse_args_cofold(sys.argv[1:])
    setup_logging(args.loglevel)
    _logger.debug("Starting reading ...")

    with open(args.cofold_out, "w") as foh:
        for line in open(args.cofold_in):
            b = line.rstrip().split(" ")
            s = b[2]
            struct = b[3]
            if args.direction == 'rgt':
                c = add_constraint(s, prime=5, struct=struct)
            else:
                c = add_constraint(s, prime=3, struct=struct)
            b[3] = c
            foh.write(" ".join(b)+"\n")

    _logger.info("Ends.")


if __name__ == '__main__':
    fix_cofold()
