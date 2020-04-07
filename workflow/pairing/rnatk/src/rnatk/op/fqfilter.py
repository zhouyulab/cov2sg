import os
import sys
import argparse
from rnatk.ios.fastq import reader_fastq
import logging
_logger = logging.getLogger(__name__)
from rnatk import __version__


def fqfilter(infile, outfile=None, minlen=1):
    assert os.path.exists(infile)
    assert minlen > 0
    outhandle = sys.stdout
    if outfile:
        outhandle = open(outfile, 'w')
    
    for fq in reader_fastq(infile):
        lenseq = len(fq.get_seq())
        lenqual = len(fq.get_qual())
        if lenseq == lenqual and lenseq >= minlen: 
            outhandle.write(str(fq)+'\n')
    
    if outfile:
        outhandle.close()


def parse_args(args):
    """Parse command line parameters"""
    parser = argparse.ArgumentParser(
        description="Parse parameters for rnatk.op.fqfilter")
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
    _logger.debug("Starting reading Fastq ...")

    fqfilter(args.fqfile, outfile=args.outfile, minlen=args.minlen)
    
    _logger.info("Ends.")


def run():
    """Entry point for console_scripts"""
    main(sys.argv[1:])
