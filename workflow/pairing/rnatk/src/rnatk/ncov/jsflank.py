#!/usr/bin/env python
"""Get JncSite flanking sequence"""
import os
from Bio import SeqIO
import sys
import argparse
import logging
_logger = logging.getLogger(__name__)


def parse_args(args):
    """Parse command line parameters"""
    parser = argparse.ArgumentParser(
        description="Parse parameters for %s" % __file__)
    parser.add_argument(
        "js",
        type=str,
        help="File for junctions sites",
        metavar="f_js")
    parser.add_argument(
        "-f",
        "--format",
        type=str,
        dest="format",
        choices=['js', 'jd'],
        help="File format: core js (js) or js by group (jd)",
        default="js")
    parser.add_argument(
        "genomefa",
        type=str,
        help="Genome Fasta file",
        metavar="genome.fa")
    parser.add_argument(
        "outfile",
        type=str,
        help="Sequences flanking JncSite",
        metavar="js.csv")
    parser.add_argument(
        "-t",
        "--top",
        type=int,
        dest="top",
        help="Top number of junctions")
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
        default="+")
    parser.add_argument(
        "-n",
        "--numevi",
        type=int,
        dest="numevi")
    parser.add_argument(
        "-l",
        "--length",
        type=int,
        dest="length",
        help="Length flanking junction site",
        default=20)
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


def parse_js(jsbed):
    for line in open(jsbed):
        b = line.rstrip().split("\t")
        chrom, start, end, name, _, strand = b[0:6]
        start = int(start)
        end = int(end)
        numevi = int(b[12])
        ntot = int(b[13])
        yield chrom, start, end, name, strand, numevi, ntot


def parse_jd(jdcsv):
    for line in open(jdcsv):
        r = line.rstrip().split("\t")
        gid, sid = r[:2]
        b = r[2:]
        chrom, start, end, name, _, strand = b[0:6]
        start = int(start)
        end = int(end)
        numevi = int(b[12])
        ntot = int(b[13])
        name = "%s|%s" % (gid, sid)
        yield chrom, start, end, name, strand, numevi, ntot


def filter_js(g_js, numevi=20, mingap=500, strand='.', top=None):
    i = 0
    for chrom, start, end, name, c_strand, c_numevi, ntot in g_js:
        if numevi is not None:
            if c_numevi < numevi:
                continue

        if end - start < mingap:
            continue

        if strand != '.':
            if c_strand != strand:
                continue

        i += 1
        if top is not None and i >= top:
            break
        yield chrom, start, end, name, c_strand, c_numevi, ntot


def g_fls(args):
    g = [r for r in SeqIO.parse(args.genomefa, "fasta")][0]
    fslen = args.length
    if args.format == "js":
        g_js = parse_js(args.js)
    elif args.format == "jd":
        g_js = parse_jd(args.js)
    for b in filter_js(g_js, numevi=args.numevi, mingap=args.mingap, strand=args.strand, top=args.top):
        chrom, start, end, name, strand, numevi, ntot = b
        fs1 = g.seq[start:start+fslen]  # rgt
        fs2 = g.seq[end:end+fslen]
        if args.direction == 'lft':
            fs1 = g.seq[start-fslen:start]
            fs2 = g.seq[end-fslen:end]

        if strand == '-':
            fs1 = fs1.reverse_complement()
            fs2 = fs2.reverse_complement()
        fs2rc = fs2.reverse_complement()
        yield chrom, start, end, name, 0, strand, numevi, ntot, fs1, fs2, fs2rc


def main(args):
    """Main entry point allowing external calls"""
    args = parse_args(args)
    setup_logging(args.loglevel)
    _logger.debug("Starting reading ...")

    with open(args.outfile, "w") as foh:
        for fls in g_fls(args):
            foh.write("\t".join(map(str, fls))+"\n")

    _logger.info("Ends.")


def run():
    """Entry point for console_scripts"""
    main(sys.argv[1:])


if __name__ == '__main__':
    run()





