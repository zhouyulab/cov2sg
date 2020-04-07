"""
Module for nCoV subgenome
"""
import argparse
import copy
import logging
import os
import sys

import pysam
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from plastid import GenomicSegment, SegmentChain
from rnatk import __version__

_logger = logging.getLogger(__name__)


class SegRead:
    """Mapped read"""
    LIBTYPE = None

    def __init__(self, r):
        assert isinstance(r, pysam.AlignedSegment)
        self.r = r

    def infer_strand(self, r):
        if SegRead.LIBTYPE is None:
            return '-' if self.r.is_reverse else '+'
        elif SegRead.LIBTYPE == "1++,1–,2+-,2-+":
            if r.is_read1:
                return '-' if self.r.is_reverse else '+'
            else:
                return '+' if self.r.is_reverse else '-'
        elif SegRead.LIBTYPE == "1+-,1-+,2++,2–":
            if r.is_read1:
                return '+' if self.r.is_reverse else '-'
            else:
                return '-' if self.r.is_reverse else '+'
        elif SegRead.LIBTYPE == "++,–":
            return '-' if self.r.is_reverse else '+'
        elif SegRead.LIBTYPE == "+-,-+":
            return '+' if self.r.is_reverse else '-'
        else:
            return "."

    @property
    def is_reverse(self):
        return self.r.is_reverse

    @property
    def name(self):
        return self.r.qname

    @property
    def chrom(self):
        return self.r.reference_name

    @property
    def strand(self):
        return self.infer_strand(self.r)

    @property
    def blocks(self):
        return self.r.blocks

    @property
    def tags(self):
        return {tag: val for tag, val in self.r.tags}

    @property
    def sgc(self):
        sgc = SegmentChain()
        for (start, end) in self.blocks:
            sgc.add_segments(GenomicSegment(self.chrom, start, end, self.strand))
        return sgc

    def __str__(self):
        return "%s:%s" % (self.name, self.sgc)

    @property
    def numjnc(self):
        return len(self.sgc.segments) - 1

    def get_jnc(self, idx):
        assert 1 <= idx <= self.numjnc
        segs = self.sgc.segments
        s = SegmentChain(*segs[:idx]).spanning_segment.end
        e = SegmentChain(*segs[idx:]).spanning_segment.start
        gs = GenomicSegment(self.chrom, s, e, self.strand)
        return Jnc(gs, self, idx)

    def g_jnc(self):
        """Generator of junctions"""
        for i in range(self.numjnc):
            yield self.get_jnc(i+1)

    def ref_loc2query_loc(self, start, end):
        """Convert reference location to corresponding location in query"""
        read = self.r
        tmp_ref_start = copy.deepcopy(read.reference_start)
        tmp_query_start = 0
        query_start = None
        query_end = None

        for cigar, length in read.cigar:
            ref_block_end = tmp_ref_start + length
            if cigar == 0:
                if tmp_ref_start <= start <= ref_block_end:
                    query_start = tmp_query_start + start - tmp_ref_start
                    query_end = query_start
            if (start is None) and (tmp_ref_start > start):
                return None

            if query_start is not None:
                if cigar == 0:
                    if ref_block_end <= end:
                        if tmp_ref_start < start < ref_block_end:
                            query_end += ref_block_end - start
                        else:
                            query_end += length
                    else:
                        query_end += max(end - max(start, tmp_ref_start), 0)
                        return query_start, query_end
                elif cigar in [1, 4, 5]:
                    query_end += length

            if cigar in [0, 1, 4, 5]:
                tmp_query_start += length
            if cigar in [0, 2, 3, 6]:
                tmp_ref_start += length

    def fetch_sequence(self, ref_start, ref_end):
        qloc = self.ref_loc2query_loc(ref_start, ref_end)
        if qloc is not None:
            query_start, query_end = qloc
            return Seq(self.r.query_sequence[query_start:query_end].upper(), IUPAC.unambiguous_dna)


class Jnc:
    """Junction from one read"""
    ARM_MINLEN = 8

    def __init__(self, gs, segread, idx):
        self.gs = gs
        self.segread = segread
        self.idx = idx

    def __str__(self):
        return "%s|%s|%d" % (str(self.gs), str(self.segread), self.idx)

    @property
    def lft(self):
        return SegmentChain(*self.segread.sgc.segments[:self.idx])

    @property
    def rgt(self):
        return SegmentChain(*self.segread.sgc.segments[self.idx:])

    @property
    def strand(self):
        return self.gs.strand

    def has_minlen(self):
        s5 = self.lft.segments[-1]
        s3 = self.rgt.segments[0]
        return min(s5.end - s5.start, s3.end - s3.start) >= Jnc.ARM_MINLEN

    @property
    def length(self):
        return self.gs.end - self.gs.start

    def shift(self, read, genome):
        """
        mapping1: AAAAA..............abcYXXBBBBB
        mapping2: AAAAAabc..............YXXBBBBB
        Convert mapping2 to mapping1
        """
        left_seq = read.fetch_sequence(self.gs.start-self.ARM_MINLEN, self.gs.start)
        right_ref_seq = genome.fetch(self.gs.chrom, self.gs.end-self.ARM_MINLEN, self.gs.end).upper()
        rev_left_seq = left_seq[::-1]
        rev_right_ref_seq = right_ref_seq[::-1]
        shift_len = 0
        for indx in range(self.ARM_MINLEN):
            if rev_left_seq[indx] == rev_right_ref_seq[indx]:
                shift_len += 1
            else:
                break
        new_gs = GenomicSegment(self.gs.chrom, self.gs.start - shift_len, self.gs.end - shift_len, self.gs.strand)
        self.gs = new_gs


class JncSite:
    """Junction Site with multiple Jnc evidences"""
    def __init__(self, jncli):
        assert len(jncli) > 0
        self.jncli = tuple(jncli)
        self.gs = jncli[0].gs
        self.don = GenomicSegment(self.gs.chrom, self.gs.start, self.gs.start+1, self.gs.strand)
        self.acc = GenomicSegment(self.gs.chrom, self.gs.end-1, self.gs.end, self.gs.strand)

    def __str__(self):
        return "\t".join(map(str, [
            self.gs.chrom, self.gs.start, self.gs.end, self.gs.strand,
            self.ntot, self.nevi]))

    @property
    def length(self):
        return self.gs.end - self.gs.start

    @property
    def ntot(self):
        return len(self.jncli)

    @property
    def nevi(self):
        """Number of evidences with non-duplicated reads"""
        return len(set(["%s+%s" % (jnc.lft, jnc.rgt) for jnc in self.jncli]))

    def dist2don(self, o):
        return max(0, max(self.don.start, o.don.start) - min(self.don.end, o.don.end))

    def dist2acc(self, o):
        return max(0, max(self.acc.start, o.acc.start) - min(self.acc.end, o.acc.end))

    def as_bed(self):
        sc = SegmentChain(self.don, self.acc)
        ts = sc.spanning_segment.start
        te = sc.spanning_segment.end
        b = sc.as_bed(thickstart=ts, thickend=te).strip("\n").split("\t")
        b.extend(map(str, [self.nevi, self.ntot]))
        return b


class JncGroup:
    """Junction Group"""
    def __init__(self, js):
        self.jsli = [js]
        self.corejs = js
        self.sc = SegmentChain(js.don, js.acc)
        self.name = None

    def distance2(self, js):
        return max(self.corejs.dist2don(js), self.corejs.dist2acc(js))

    def add(self, js):
        self.jsli.append(js)

    def merge(self):
        for i in range(1, len(self.jsli)):
            js = self.jsli[i]
            self.sc.add_segments(js.don)
            self.sc.add_segments(js.acc)

    def get_ntot(self):
        return sum(js.ntot for js in self.jsli)

    def get_nevi(self):
        return sum(js.nevi for js in self.jsli)

    def set_name(self, name):
        self.name = name

    def as_bed(self):
        if self.name is not None:
            self.sc.attr['ID'] = self.name
        if self.sc.strand == "+":
            color = "255,0,0"
        else:
            color = "0,0,255"
        self.sc.attr['color'] = color
        b = self.sc.as_bed().strip("\n").split("\t")
        b.extend(map(str, [self.get_nevi(), self.get_ntot()]))
        return b

    def __str__(self):
        return "\t".join(map(str, [
            self.corejs, self.get_ntot(), self.get_nevi()]))

    def jncsites(self):
        jncs = []
        for i in range(len(self.jsli)):
            jncs.append([self.name, str(i+1)] + self.jsli[i].as_bed())
            i += 1
        return jncs


class Genome:
    """
    Small genome as dict of chrom:Seq
    """
    def __init__(self, fagenome):
        assert os.path.exists(fagenome)
        self.fa = fagenome
        self.c2s = self._read_genome()

    def _read_genome(self):
        c2s = {}
        for r in SeqIO.parse(self.fa, "fasta"):
            c2s[r.id] = r.seq.upper()
        return c2s

    def fetch(self, chrom, start, end):
        if chrom not in self.c2s:
            return None
        return self.c2s[chrom][start:end]


class JncAssembler:
    """
    Junction assembler
    """
    JID = 0

    def __init__(self, bam, window=5, minevi=2, mingap=30, genome=None):
        self.bam = bam
        self.window = window
        self.minevi = minevi
        self.mingap = mingap
        self.gd = genome
        if genome is not None:
            assert os.path.exists(genome)
            self.gd = Genome(genome)
        self.jsdi = {}
        self.jgli = []

    def numjs(self):
        return len(self.jsdi)

    def numjg(self):
        return len(self.jgli)

    def sorted(self):
        for gs, js in sorted(self.jsdi.items(), key=lambda x: (-x[1].nevi, -x[1].ntot)):
            yield gs, js

    def g_segread(self, chrom, start, end):
        for r in self.bam.fetch(chrom, start, end):
            if r.is_secondary:
                continue
            yield SegRead(r)

    def has_arm_match(self, read, chrom, js_start_site, js_end_site, arm_len):
        """Require both arms not having mismatches"""
        read_left_seq = read.fetch_sequence(js_start_site - arm_len, js_start_site)
        if read_left_seq is None:  # skip if arm is not long enough
            return False
        read_right_seq = read.fetch_sequence(js_end_site, js_end_site + arm_len)
        if read_right_seq is None:
            return False
        ref_left_seq = self.gd.fetch(chrom, js_start_site - arm_len, js_start_site)
        ref_right_seq = self.gd.fetch(chrom, js_end_site, js_end_site + arm_len)
        if read.is_reverse:
            read_left_seq = read_left_seq.reverse_complement()
            read_right_seq = read_right_seq.reverse_complement()
        return read_left_seq == ref_left_seq and read_right_seq == ref_right_seq

    def get_jncs(self, chrom, start, end, strand):
        d = {}
        for sr in self.g_segread(chrom, start, end):
            if sr.strand != strand:
                continue

            for jnc in sr.g_jnc():
                if not (start <= jnc.gs.start and jnc.gs.end <= end):
                    continue
                if not jnc.has_minlen():
                    continue
                if self.gd is not None:
                    if not self.has_arm_match(sr, chrom, jnc.gs.start, jnc.gs.end, Jnc.ARM_MINLEN):
                        continue

                if jnc.gs not in d:
                    d[jnc.gs] = [jnc]
                else:
                    d[jnc.gs].append(jnc)

        self.jsdi = {}
        for gs, jncli in d.items():
            self.jsdi[gs] = JncSite(jncli)

    def cluster(self):
        self.jgli = []
        for gs, js in self.sorted():  # greedy from strongest junctions
            if js.length < self.mingap:  # filter out small junctions
                continue

            if js.nevi < self.minevi:  # filter out junctions with low support
                continue

            found_jg = False
            for jg in self.jgli:
                dist = jg.distance2(js)
                if dist <= self.window:
                    found_jg = True
                    jg.add(js)
                    break

            if not found_jg:
                new_jg = JncGroup(js)
                self.jgli.append(new_jg)

    def merge(self):
        for jg in self.jgli:
            jg.merge()

    def assemble(self, chrom, start, end, strand):
        _logger.info("\nFinding junctions ...")
        self.get_jncs(chrom, start, end, strand)

        _logger.info("\nNumber JncSite: %d", self.numjs())
        if self.numjs() == 0:
            return []

        _logger.info("\nClustering junctions to groups ...")
        self.cluster()
        _logger.info("\nNumber JncGroup: %d", self.numjg())

        _logger.info("\nMerging in JncGroup ...")
        self.merge()

    def g_jnc_beds(self):
        for jg in self.jgli:
            JncAssembler.JID += 1
            jg.set_name(str(JncAssembler.JID))
            jgbed = jg.as_bed()
            jsbed = jg.corejs.as_bed()
            for ix in (3, 8):  # set name
                jsbed[ix] = jgbed[ix]
            yield jgbed, jsbed, jg.jncsites()

    def i_chromsize(self):
        for chrom, size in zip(self.bam.references, self.bam.lengths):
            yield chrom, size

    def assemble_genome(self, outprefix):
        with open(outprefix + "_js.bed", "w") as foh_js:
            with open(outprefix + "_jg.bed", "w") as foh_jg:
                with open(outprefix + "_jd.csv", "w") as foh_jd:
                    for chrom, size in self.i_chromsize():
                        start = 0
                        for strand in ("+", "-"):
                            _logger.info("\nAnalzye %s:%d-%d:%s ..." % (chrom, start, size, strand))
                            self.assemble(chrom, start, size, strand)
                            for jgbed, jsbed, jncsites in self.g_jnc_beds():
                                foh_jg.write("\t".join(map(str, jgbed)) + "\n")
                                foh_js.write("\t".join(map(str, jsbed)) + "\n")
                                for jncsite in jncsites:
                                    foh_jd.write("\t".join(map(str, jncsite)) + "\n")


def parse_args(args):
    """Parse command line parameters"""
    parser = argparse.ArgumentParser(
        description="Parse parameters for rnatk.ncov.subgenome")
    parser.add_argument(
        "fbam",
        type=str,
        help="BAM file with index",
        metavar="infile.bam")
    parser.add_argument(
        "--genome",
        type=str,
        dest="genome",
        help="genome.fa")
    parser.add_argument(
        "-w",
        "--window",
        type=int,
        dest="window",
        help="Window size to merge nearby junctions",
        default=5)
    parser.add_argument(
        "-m",
        "--minevi",
        type=int,
        dest="minevi",
        help="Minimum number of non-duplicated reads evidence",
        default=1)
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
        "-r",
        "--rule",
        type=str,
        dest="rule",
        help="Library rule to infer read strand",
        default="1+-,1-+,2++,2–")
    parser.add_argument(
        "-o",
        "--outprefix",
        type=str,
        help="Output prefix for _js.bed, _jg.bed")
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
    _logger.debug("Starting assembling junctions ...")
    SegRead.LIBTYPE = args.rule
    if args.armsize is not None:
        Jnc.ARM_MINLEN = args.armsize
    bam = pysam.AlignmentFile(args.fbam, "rb")
    ja = JncAssembler(bam, window=args.window, minevi=args.minevi, mingap=args.mingap, genome=args.genome)
    ja.assemble_genome(args.outprefix)

    _logger.info("Ends.")


def run():
    """Entry point for console_scripts"""
    main(sys.argv[1:])


if __name__ == "__main__":
    run()
