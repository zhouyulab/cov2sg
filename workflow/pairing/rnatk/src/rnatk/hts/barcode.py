"""
Extract barcode/index from read and rename its name as num|barcode
"""

import os
import sys
import argparse
from rnatk.ios.fastq import reader_fastq, Fastq
import logging
_logger = logging.getLogger(__name__)
from rnatk import __version__


def iter_parsed_read(fqfile, bc_start=0, bc_end=4, bc_drop=True):
    """Generator of parsed Fastq record"""
    i = 0
    for fq in reader_fastq(fqfile):
        i += 1
        old_seq = fq.get_seq()
        old_qual = fq.get_qual()
        barcode = old_seq[bc_start:bc_end].upper()
        
        new_name = "%d|%s" % (i, barcode) # rename read
        new_seq = old_seq
        new_qual = old_qual
        if bc_drop:
            new_seq = old_seq[:bc_start] + old_seq[bc_end:]
            new_qual = old_qual[:bc_start] + old_qual[bc_end:]
            
        yield Fastq(new_name, new_seq, new_qual)


def parse_barcode(args):
    foh = sys.stdout
    if args.outfile:
        foh = open(args.outfile, "w")

    for fq in iter_parsed_read(
        args.fqfile, bc_start=args.bc_start, bc_end=args.bc_end, bc_drop=args.drop):
        if fq.get_length() >= args.minlen:
            foh.write(str(fq)+"\n")

    if args.outfile:
        foh.close()


def parse_args(args):
    """Parse command line parameters"""
    parser = argparse.ArgumentParser(
        description="Parse parameters for rnatk.parse_barcode")
    parser.add_argument(
        "fqfile",
        type=str,        
        help="Fastq file",
        metavar="infile.fq")
    parser.add_argument(
        "-o",
        "--outfile", 
        type=str, 
        dest="outfile",
        help="Outfile name or sys.stdout")
    parser.add_argument(
        "-s",
        "--start",
        type=int,        
        dest="bc_start",
        help="0-based start position",
        default=0)
    parser.add_argument(
        "-e",
        "--end",
        type=int,        
        dest="bc_end",
        help="1-based end position",
        default=4)
    parser.add_argument(
        "--drop",
        action="store_true",
        help="Drop barcode sequence",
        )
    parser.add_argument(
        "-m",
        "--minlen",
        type=int,        
        dest="minlen",
        help="Minimum length of parsed read",
        default=1)
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
    _logger.debug("Starting parsing barcode ...")

    parse_barcode(args)

    _logger.info("Ends.")


def run():
    """Entry point for console_scripts"""
    main(sys.argv[1:])



