#!/usr/bin/env python
import os
import sys
import logging
_logger = logging.getLogger(__name__)
import argparse
from rnatk import __version__
from rnatk.ios.fastq import reader_fastq
from rnatk.util.rc import reverse_complement, reverse
from Bio import pairwise2


class PEMerger(object):
    """
    Pair-End reader correction and merger
    """
    def __init__(self, r1fq, r2fq, outprefix, rawlen=150, min_match_ratio=0.7):
        assert os.path.exists(r1fq)
        assert os.path.exists(r2fq)
        self.r1fq = r1fq
        self.r2fq = r2fq
        self.outprefix = outprefix
        self.rawlen = rawlen
        self.min_match_ratio = min_match_ratio

    def pemerge(self):
        """
        Merge Pair-end reads by semi-global overlapping alignment
        """
        for r1, r2 in zip(reader_fastq(self.r1fq), reader_fastq(self.r2fq)):
            s1 = r1.get_seq()
            s2 = r2.get_seq()
            q1 = r1.get_qual()
            q2 = r2.get_qual()
            s2rc = reverse_complement(s2)
            q2rc = reverse(r2.get_qual())
            s1c, q1c, s2rcc, q2rcc, match_ratio = align_and_correct(
                s1, q1, s2rc, q2rc, min_match_ratio=self.min_match_ratio)
            yield ((r1.get_name(), s1c, q1c),
                   (r2.get_name(), reverse_complement(s2rcc), reverse(q2rcc)),
                   match_ratio)

    def merge(self):
        with open(self.outprefix+"_R1.fq", "w") as foh1, open(self.outprefix+"_R2.fq", "w") as foh2:
            for s1, s2, match_ratio in self.pemerge():
                foh1.write("\n".join(["@%s" % s1[0], s1[1], "+", s1[2]]) + "\n")
                foh2.write("\n".join(["@%s" % s2[0], s2[1], "+", s2[2]]) + "\n")


def align_and_correct(s1, q1, s2, q2, min_match_ratio=0.7):
    """
    Align and correct two sequences with qualities: (s1, q1) and (s2, q2)
    """
    alignment = semiglobal(s1, s2)
    (aln1, aln2, score, begin, end) = alignment
    assert begin == 0 and end >= len(s1) and end >= len(s2)
    nmatch = 0
    for i in range(0, len(aln1)):
        if aln1[i] == aln2[i]:
            nmatch += 1
    match_ratio = nmatch / min(len(s1), len(s2))
    if match_ratio < min_match_ratio: # skip correction 
        return s1, q1, s2, q2, match_ratio

    cor1, cor2 = [], []
    idx1, idx2 = -1, -1
    for i in range(0, len(aln1)):
        if aln1[i] != "-":
            idx1 += 1
        if aln2[i] != "-":
            idx2 += 1
        if aln1[i] != "-" and aln2[i] != "-":
            if aln1[i] != aln2[i]:
                if q1[idx1] >= q2[idx2]:
                    b = (s1[idx1], q1[idx1])
                else:
                    b = (s2[idx2], q2[idx2])
            else:
                b = (s1[idx1], max(q1[idx1], q2[idx2]))
            cor1.append(b)
            cor2.append(b)
        elif aln1[i] != "-" and aln2[i] == "-":
            cor1.append((s1[idx1], q1[idx1]))
        elif aln1[i] == "-" and aln2[i] != "-":
            cor2.append((s2[idx2], q2[idx2]))
    
    s1c = "".join([b[0] for b in cor1])
    q1c = "".join([b[1] for b in cor1])
    s2c = "".join([b[0] for b in cor2])
    q2c = "".join([b[1] for b in cor2])
    return s1c, q1c, s2c, q2c, match_ratio


def semiglobal(s1, s2):
    """Do semi-global alignment"""
    alignments = pairwise2.align.globalms(
        s1, s2, 1, -1, -3, -1,
        penalize_end_gaps=(False, False),
        one_alignment_only=True)
    return alignments[0]


def parse_args(args):
    """Parse command line parameters"""
    parser = argparse.ArgumentParser(
        description="Parse parameters for rnatk.op.pemerge")
    parser.add_argument(
        "r1fq",
        type=str,        
        help="R1 fastq file",
        metavar="name_R1.fq")
    parser.add_argument(
        "r2fq",
        type=str,        
        help="R2 fastq file",
        metavar="name_R2.fq")
    parser.add_argument(
        "-o",
        "--outprefix", 
        type=str, 
        dest="outprefix",
        help="Outprefix of multiple output files")
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
    _logger.debug("Starting processing ...")

    pm = PEMerger(args.r1fq, args.r2fq, outprefix=args.outprefix)
    pm.merge()
    
    _logger.info("Ends.")


def run():
    """Entry point for console_scripts"""
    main(sys.argv[1:])

