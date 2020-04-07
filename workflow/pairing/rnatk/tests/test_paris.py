#!/usr/bin/env python
"""Tester for paris module"""
import pysam
import pytest
from plastid import GenomicSegment, SegmentChain
from rnatk.rna.paris import (ParisRead, find_duplex, sgc2str, sgc_distance,
                             sgc_overlap, sgc_union)


@pytest.fixture
def sgc1():
    exon1 = GenomicSegment("chr1", 100, 200, "+")
    exon2 = GenomicSegment("chr1", 300, 400, "+")
    return SegmentChain(exon1, exon2, ID="SGC1", alias="SGC1")


@pytest.fixture
def sgc2():
    exon1 = GenomicSegment("chr1", 1150, 1200, "+")
    exon2 = GenomicSegment("chr1", 1300, 1450, "+")
    return SegmentChain(exon1, exon2, ID="SGC2", alias="SGC2")


@pytest.fixture
def sgc3():
    return SegmentChain(GenomicSegment("chr1", 350, 450, "+"))


def test_segmentchain(sgc1):
    sgc1.add_segments(GenomicSegment("chr1", 350, 450, "+"))
    assert str(sgc1) == "chr1:100-200^300-450(+)"


def test_sgc2str(sgc1):
    assert sgc2str(sgc1) == "chr1:100-200^300-400(+)[2]"


def test_sgc_overlap_no(sgc1, sgc2):
    newsc = sgc_overlap(sgc1, sgc2)
    assert newsc.length == 0 and newsc.segments == []


def test_sgc_overlap_yes(sgc1, sgc3):
    newsc = sgc_overlap(sgc1, sgc3)
    assert str(newsc) == "chr1:350-400(+)"


def test_sgc_union(sgc1, sgc2):
    newsc = sgc_union(sgc1, sgc2)
    assert str(newsc) == "chr1:100-200^300-400^1150-1200^1300-1450(+)"


def test_sgc_overlap_yes(sgc1, sgc3):
    newsc = sgc_union(sgc1, sgc3)
    assert str(newsc) == "chr1:100-200^300-450(+)"


def test_sgc_distance_nooverlap(sgc1, sgc2):
    d = sgc_distance(sgc1, sgc2)
    assert d == 750


def test_sgc_distance_overlap(sgc1, sgc3):
    d = sgc_distance(sgc1, sgc3)
    assert d == 0


def test_find_duplex():
    bam = pysam.AlignmentFile("tests/data/trsread.bam", "rb")
    rli = find_duplex(bam, "MN996528", 0, 29891)
    assert len(rli) == 10


class TestParisRead:
    def test_init_and_gduplex(self):
        bam = pysam.AlignmentFile("tests/data/trsread.bam", "rb")
        rli = find_duplex(bam, "MN996528", 0, 29891)
        r = [r for r in rli if r.qname == '1'][0]
        pr = ParisRead(r)
        assert str(pr) == "1\tMN996528:105-130^330-355(+)"
        dpli = list(pr.g_duplex())
        assert len(dpli) == 1
        assert str(dpli[0]) == "MN996528:105-130(+)[1] & MN996528:330-355(+)[1]"
