import os
import pytest
import pysam
from rnatk.bin.csv2filter import csv2filter


@pytest.fixture()
def outcsv():
    tmpcsv = "tests/data/csv2filter.tmp.csv"
    yield tmpcsv
    os.unlink(tmpcsv)


def test_csv2filter_num1(outcsv):
    csv2filter("tests/data/dupbc.bed", '<', 0, col=5, outfile=outcsv)
    lines = [line for line in open(outcsv)]
    assert len(lines) == 0


def test_csv2filter_num2(outcsv):
    csv2filter("tests/data/dupbc.bed", '>=', 0, col=5, outfile=outcsv)
    lines = [line for line in open(outcsv)]
    assert len(lines) == 16


def test_csv2filter_num3(outcsv):
    csv2filter("tests/data/dupbc.bed", '!=', 0, col=5, outfile=outcsv)
    lines = [line for line in open(outcsv)]
    assert len(lines) == 0


def test_csv2filter_str1(outcsv):
    csv2filter("tests/data/dupbc.bed", '==', 'chr2', col=1, outfile=outcsv)
    lines = [line for line in open(outcsv)]
    assert len(lines) == 2


def test_csv2filter_err1(outcsv):
    isok = csv2filter("tests/data/dupbc.bed", '==', 5, col=1, outfile=outcsv)
    lines = [line for line in open(outcsv)]
    assert len(lines) == 0 and isok is False
