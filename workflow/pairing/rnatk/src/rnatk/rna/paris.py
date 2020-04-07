#!/usr/bin/env python
"""Module for PARIS data"""
import math
import copy
import collections
import pysam
from plastid import GenomicSegment, SegmentChain
from operator import attrgetter


def is_overlap(iv1, iv2):
    """Check if one interval (start, end) overlaps with the other"""
    return min(iv1[1], iv2[1]) - max(iv1[0], iv2[0]) > 0


def find_duplex(bam, chrom, start, end):
    """Find duplex read whose segments are in given window"""
    rli = []
    for r in bam.fetch(chrom, start, end):
        isin = True
        for b in r.get_blocks():
            if not is_overlap((start, end), b):
                isin = False
        if isin:
            rli.append(r)

    return rli


def duplex2bam(bam, chrom, start, end, fout):
    """Write duplex within given window into BAM file"""
    outbam = pysam.AlignmentFile(fout, "wb", template=bam)
    rli = find_duplex(bam, chrom, start, end)
    for r in rli:
        outbam.write(r)

    outbam.close()
    return len(rli)


def sgc2str(sgc):
    """Append number of segments to SegmentChain"""
    return "%s[%d]" % (sgc, len(sgc.segments))


def sgc_overlap(sgc1, sgc2):
    """Get overlapped SegmentChain between given two"""
    sgc1.add_masks(*sgc2.segments)
    new_sgc = sgc1.get_masks_as_segmentchain()
    sgc1.reset_masks()
    return new_sgc


def sgc_union(sgc1, sgc2):
    new_sgc = copy.deepcopy(sgc1)
    new_sgc.add_segments(*sgc2.segments)
    return new_sgc


def sgc_distance(sgc1, sgc2):
    s1 = sgc1.spanning_segment.start
    s2 = sgc2.spanning_segment.start
    e1 = sgc1.spanning_segment.end
    e2 = sgc2.spanning_segment.end
    return max(0, max(s1, s2) - min(e1, e2))


class ParisRead:
    """Read from PARIS-seq"""
    def __init__(self, r):
        assert isinstance(r, pysam.AlignedSegment)
        self.r = r
        self.chrom = r.reference_name
        self.strand = '-' if self.r.is_reverse else '+'
        self.blocks = self.r.blocks
        self.tags = {tag: val for tag, val in self.r.tags}
        self.sgc = SegmentChain()
        for (start, end) in self.r.blocks:
            self.sgc.add_segments(
                GenomicSegment(self.chrom, start, end, self.strand))

    def __str__(self):
        return "%s\t%s" % (self.r.qname, self.sgc)

    def g_duplex(self):
        """Generator of duplex by gaps"""
        segs = self.sgc.segments
        for i in range(1, len(segs)):
            lftsgc = SegmentChain(*segs[:i], ID="%dLeft" % i)
            rgtsgc = SegmentChain(*segs[i:], ID="%dRight" % i)
            yield Duplex(lftsgc, rgtsgc, self)


class Duplex:
    """Duplex with left and right arms from PARIS read"""
    ARM_MINLEN = 5

    def __init__(self, lsgc, rsgc, pr=None):
        self.larm = lsgc
        self.rarm = rsgc
        self.pr = pr
        self.start = self.larm.spanning_segment.start
        self.end = self.rarm.spanning_segment.end

    def __str__(self):
        return "%s & %s" % (sgc2str(self.larm), sgc2str(self.rarm))

    def is_spliced(self, sjdb=None):
        if sjdb is not None:
            jc = self.get_junc()
            return sjdb.has_splicejunc(jc)

        return False

    def get_junc(self):
        js = self.larm.spanning_segment.end
        je = self.rarm.spanning_segment.start
        return GenomicSegment(self.larm.chrom, js, je, self.larm.strand)

    def has_minlen(self):
        return min(self.lenlarm(), self.lenrarm()) >= Duplex.ARM_MINLEN

    def intersect(self, other):
        larm = sgc_overlap(self.larm, other.larm)
        rarm = sgc_overlap(self.rarm, other.rarm)
        return Duplex(larm, rarm)

    def union(self, other):
        larm = sgc_union(self.larm, other.larm)
        rarm = sgc_union(self.rarm, other.rarm)
        return Duplex(larm, rarm)

    def distance(self, other):
        return (sgc_distance(self.larm, other.larm),
                sgc_distance(self.rarm, other.rarm))

    def lenlarm(self):
        return self.larm.length

    def lenrarm(self):
        return self.rarm.length

    def mergable(self, other):
        dp = self.union(other)
        if (max(self.distance(other)) < 10 and
                max(dp.larm.length, dp.rarm.length) < 30):
            return True
        return False

    def merge(self, other):
        return self.union(other)


class DuplexGroup:
    """Duplex Group"""

    def __init__(self, dpx):
        self.dpli = [dpx]
        self.core = dpx
        self.conn_score = None
        self.id = None

    def __str__(self):
        dgstr = "%s=%d" % (self.core, len(self.dpli))
        if self.conn_score is not None:
            dgstr += "=%.4f" % self.conn_score
        if self.id is not None:
            dgstr += "=%s" % str(self.id)
        return dgstr

    def add(self, dpx, new_core):
        self.dpli.append(dpx)
        self.core = new_core

    def intersect_duplex(self, dpx):
        return self.core.intersect(dpx)

    def mergable(self, other):
        return self.core.mergable(other.core)

    def merge(self, other):
        self.dpli += other.dpli
        self.core = self.core.merge(other.core)

    def num_unique(self):
        return len(set([
            (dp.start, dp.end) for dp in self.dpli]))

    def num_dp(self):
        return len(self.dpli)

    def set_score(self, score):
        self.conn_score = score

    def set_id(self, dgid):
        self.id = dgid


class DGA:
    """Duplex Group Assembler"""
    DGnum = 0

    def __init__(self, bam, sjdb):
        self.bam = bam
        self.sjdb = sjdb
        self.numr = 0
        self.numdp = 0
        self.numdprm = 0

    def novdg(self, dgli):
        """Non-overlapping DG"""
        pass

    def assemble(self, chrom, start, end, strand):
        dpli = self.get_duplex(chrom, start, end, strand)
        print("\nNumber duplex: ", len(dpli))
        if len(dpli) == 0:
            return []

        dgli = self.cluster(dpli)
        print("\nNumber duplex group: ", len(dgli))
        for dg in dgli:
            print(dg)
        if len(dgli) == 0:
            return []

        dgmerged = self.merge(dgli)
        print("\nNumber merged duplex group: ", len(dgmerged))
        for dg in dgmerged:
            print(dg)

        dgfi = [dg for dg in dgmerged if self.screen_dg(dg)]
        print("\nNumber final duplex group: ", len(dgfi))

        for dg in sorted(dgfi, key=lambda x: -x.num_dp()):
            self.add_dgid(dg)
            print(dg)
        return dgfi

    def add_dgid(self, dg):
        DGA.DGnum += 1
        dg.set_id('%06d' % DGA.DGnum)

    def screen_dg(self, dg):
        conn_score = self.score_dg(dg)
        if dg.num_unique() >= 2 and conn_score >= 0.01:
            return True
        return False

    def score_dg(self, dg):
        a = self.coverage(dg.core.larm)
        b = self.coverage(dg.core.rarm)
        conn_score = len(dg.dpli) / math.sqrt(a * b)
        dg.set_score(conn_score)
        return conn_score

    def coverage(self, sgc):
        chrom = sgc.chrom
        start = sgc.spanning_segment.start
        end = sgc.spanning_segment.end
        strand = sgc.strand
        n = 0
        for r in self.bam.fetch(chrom, start, end):
            if not r.is_unmapped:
                r_strand = '-' if r.is_reverse else '+'
                if r_strand == strand:
                    n += 1
        return n

    def merge(self, dgli):
        dgsorted = list(sorted(dgli, key=lambda s: (s.core.start, s.core.end)))
        for dg in dgli:
            assert isinstance(dg, DuplexGroup), dg
        queue = collections.deque(dgsorted)
        dgfi = []
        dg = queue.popleft()
        while len(queue) > 0:
            dg2 = queue.popleft()
            if dg.mergable(dg2):
                dg.merge(dg2)
            else:
                dgfi.append(dg)
                dg = dg2

        dgfi.append(dg)
        return dgfi

    def cluster(self, dpli):
        dpsorted = list(sorted(dpli, key=attrgetter('start', 'end')))
        dgli = []
        for dpx in dpli:
            found_dg = False
            for dg in dgli:
                overlapped = dg.intersect_duplex(dpx)
                if overlapped.has_minlen():
                    found_dg = True
                    dg.add(dpx, overlapped)
                    break

            if not found_dg:
                new_dg = DuplexGroup(dpx)
                dgli.append(new_dg)

        return dgli

    def get_duplex(self, chrom, start, end, strand):
        self.numr = 0
        self.numdp = 0
        self.numdprm = 0
        dpli = []
        for r in find_duplex(self.bam, chrom, start, end):
            r_strand = '-' if r.is_reverse else '+'
            if r_strand != strand:
                continue

            self.numr += 1
            for dp, keep in self.g_duplex_from_read(r):
                self.numdp += 1
                if keep:
                    dpli.append(dp)
                else:
                    self.numdprm += 1
        return dpli

    def g_duplex_from_read(self, r):
        pr = ParisRead(r)
        for dp in pr.g_duplex():
            keep = (not dp.is_spliced(self.sjdb)) and dp.has_minlen()
            yield dp, keep
