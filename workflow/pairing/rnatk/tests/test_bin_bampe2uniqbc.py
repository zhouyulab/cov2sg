import os
import pytest
import pysam
from rnatk.bin.bampe2uniqbc import bampe2uniqbc


@pytest.fixture()
def outbam():
    tmpbam = "tests/data/dupbcpe.expected.bam"
    yield tmpbam
    os.unlink(tmpbam)


def test_bampe2uniqbc1(outbam):
    inbam = "tests/data/dupbcpe.bam"
    bampe2uniqbc(inbam, outbam, maxtop=1, bclen=4)
    bam = pysam.Samfile(outbam, "rb")
    rns = [r.qname for r in bam.fetch(until_eof=True)]
    assert len(set(rns)) == 12


def test_bampe2uniqbc3(outbam):
    inbam = "tests/data/dupbcpe.bam"
    bampe2uniqbc(inbam, outbam, maxtop=3, bclen=4)
    bam = pysam.Samfile(outbam, "rb")
    rns = [r.qname for r in bam.fetch(until_eof=True)]
    assert len(set(rns)) == 16
