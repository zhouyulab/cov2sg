#!/usr/bin/env python
"""
Convert Bed 8 to Bed 12 format
"""
import os
import sys
import argparse
from rnatk.ucsc.track import TrackBed8
import logging
_logger = logging.getLogger(__name__)
from rnatk import __version__

def bed8To12(infile, outfile, rgb_fwd='0,0,0', rgb_rev='0,0,0'):
    tmk = TrackBed8(rgb_fwd, rgb_rev)
    tmk.bedTo12(infile, outfile, sort=False)


def parse_args(args):
    """Parse command line parameters"""
    parser = argparse.ArgumentParser(
        description="Parse parameters for rnatk.bin.bed8To12")
    parser.add_argument(
        "infile",
        type=str,        
        help="Input file",
        metavar="in.bed")
    parser.add_argument(
        "outfile",
        type=str,        
        help="Output file",
        metavar="out.bed")
    parser.add_argument(
        "--rgbfwd",
        dest="rgbfwd",
        type=str,
        default='0,0,0',
        help="RGB for Forward strand (default '0,0,0')")
    parser.add_argument(
        "--rgbrev",
        dest="rgbrev",
        type=str,
        default='0,0,0',
        help="RGB for Reverse strand (default '0,0,0')")
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
    _logger.debug("Starting reading BED ...")

    bed8To12(args.infile, args.outfile, args.rgbfwd, args.rgbrev)

    _logger.info("Ends.")

def run():
    """Entry point for console_scripts"""
    main(sys.argv[1:])

if __name__ == "__main__":
    run()


