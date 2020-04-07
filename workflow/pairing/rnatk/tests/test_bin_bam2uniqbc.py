import os
import pytest
import pysam
from rnatk.bin.bam2uniqbc import bam2uniqbc


@pytest.fixture()
def outbam():
    tmpbam = "tests/data/dupbc.expected.bam"
    yield tmpbam
    os.unlink(tmpbam)


def test_bam2uniqbc1(outbam):
    inbam = "tests/data/dupbc.bam"
    bam2uniqbc(inbam, outbam, maxtop=1, bclen=4)
    bam = pysam.Samfile(outbam, "rb")
    rns = [r.qname for r in bam.fetch(until_eof=True)]
    assert len(set(rns)) == 12


def test_bam2uniqbc3(outbam):
    inbam = "tests/data/dupbc.bam"
    bam2uniqbc(inbam, outbam, maxtop=3, bclen=4)
    bam = pysam.Samfile(outbam, "rb")
    rns = [r.qname for r in bam.fetch(until_eof=True)]
    assert len(set(rns)) == 16

