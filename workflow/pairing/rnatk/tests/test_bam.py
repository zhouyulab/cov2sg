import pytest
import pysam
from rnatk.ios.bam import *


@pytest.fixture
def bamin():
    return pysam.Samfile("tests/data/dupbc.bam", "rb")


def test_ibam_chrom(bamin):
    c_n = [(ci[0], len(list(ci[1]))) for ci in ibam_chrom(bamin)]
    assert c_n == [('chr1', 14), ('chr2', 2)]


def test_ibam_chromstart(bamin):
    c_r = [(ci[0], len(ci[1])) for ci in ibam_chromStart(bamin)]
    assert c_r == [(('chr1', 10), 14), (('chr2', 5), 2)]


def test_igroup_strand(bamin):
    c_r = []
    for chrom, rli in ibam_chromStart(bamin):
        for strand, reads in igroup_strand(rli):
            c_r.append((chrom, strand, len(reads)))
    assert c_r == [
        (('chr1', 10), '+', 7),
        (('chr1', 10), '-', 7),
        (('chr2', 5), '+', 1),
        (('chr2', 5), '-', 1)]


def test_igroup_end(bamin):
    c_r = []
    for chrom, rli in ibam_chrom(bamin):
        for end, reads in igroup_end(rli):
            c_r.append((chrom, end, len(reads)))
    assert c_r == [
        ('chr1', 25, 2),
        ('chr1', 30, 10),
        ('chr1', 35, 2),
        ('chr2', 20, 2)]


def test_ibam_chromstartendstrand(bamin):
    c_r = [(ci[0], len(ci[1])) for ci in ibam_chromStartEndStrand(bamin)]
    assert c_r == [
        (('chr1', 10, 25, '+'), 1),
        (('chr1', 10, 25, '-'), 1),
        (('chr1', 10, 30, '+'), 5),
        (('chr1', 10, 30, '-'), 5),
        (('chr1', 10, 35, '+'), 1),
        (('chr1', 10, 35, '-'), 1),
        (('chr2', 5, 20, '+'), 1),
        (('chr2', 5, 20, '-'), 1), ]


@pytest.fixture
def bampein():
    return pysam.Samfile("tests/data/dupbcpe.bam", "rb")


def test_ibampe(bampein):
    c_n = [(ci[0], len(ci[1])) for ci in ibampe(bampein)]
    assert c_n == [('chr1', 14), ('chr2', 2)]


def test_ibampe_chromstartendr1strand(bampein):
    s_rli = [ci for ci in ibampe_chromStartEndR1Strand(bampein)]
    s_rns = []
    for sig, rli in s_rli:
        s_rns.append((sig, tuple([pe['R1'].qname for pe in rli])))
    assert s_rns == [
        (('chr1', 10, 125, '+'), ('f6|AAAA',)),
        (('chr1', 10, 125, '-'), ('r6|AAAA',)),
        (('chr1', 10, 130, '+'), ('f1|AAAA', 'f2|AAAA', 'f3|AAAA', 'f4|CCCC', 'f5|GGGG')),
        (('chr1', 10, 130, '-'), ('r1|AAAA', 'r2|AAAA', 'r3|AAAA', 'r4|CCCC', 'r5|GGGG')),
        (('chr1', 10, 135, '+'), ('f7|AAAA',)),
        (('chr1', 10, 135, '-'), ('r7|AAAA',)),
        (('chr2', 5, 120, '+'), ('2f|AAAA',)),
        (('chr2', 5, 120, '-'), ('2r|AAAA',))]
