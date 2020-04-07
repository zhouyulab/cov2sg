#!/usr/bin/env python
"""
Remove reads having same location (chr, start, end, strand) by barcode sequences.

barcode/index sequence is append to the end of read name as `read|barcode`.
if barcode length is set as 0, reads are collapsed by position only.
Input BAM must be sorted by position with samtools sort.
"""
import os
import sys
import argparse
import pysam
from rnatk.ios.bam import ibam_chromStartEndStrand, igroup_barcode
import logging
_logger = logging.getLogger(__name__)
from rnatk import __version__


def bam2uniqbc(inbam, outbam, maxtop=1, bclen=4):
    assert os.path.exists(inbam)
    bamin = pysam.Samfile(inbam, "rb")
    oh = pysam.Samfile(outbam, "wb", template=bamin)
    for sig, reads in ibam_chromStartEndStrand(bamin):
        for bc, grpreads in igroup_barcode(reads, bclen):
            for s in grpreads[:maxtop]:
                oh.write(s)
    oh.close()


def parse_args(args):
    """Parse command line parameters"""
    parser = argparse.ArgumentParser(
        description="Parse parameters for rnatk.bin.bam2uniqbc")
    parser.add_argument(
        "inbam",
        type=str,        
        help="Input Bam file, position sorted with samtools",
        metavar="input.bam")
    parser.add_argument(
        "outbam",
        type=str,        
        help="Output Bam file",
        metavar="output.bam")
    parser.add_argument(
        #"-t",
        "--top",
        type=int,        
        dest="top",
        help="Maximum number of reads to keep[1]",
        default=1)
    parser.add_argument(
        #"-l",
        "--bclen",
        type=int,        
        dest="bclen",
        help="Length of barcode[4]",
        default=4)
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
    _logger.debug("Starting reading Bam ...")

    bam2uniqbc(args.inbam, args.outbam, maxtop=1, bclen=4)
    
    _logger.info("Ends.")

def run():
    """Entry point for console_scripts"""
    main(sys.argv[1:])

if __name__ == "__main__":
    run()

