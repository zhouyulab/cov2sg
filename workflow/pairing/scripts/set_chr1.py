#!/usr/bin/env python
"""
Reset chrom as chr1 for visualization in UCSC genome browser
"""
import os
import sys


if len(sys.argv) < 2:
    sys.stderr.write("Usage: program in.bed > out.bed\n")
    sys.exit(1)

fbed = sys.argv[1]
assert os.path.exists(fbed), fbed
for line in open(fbed):
    b = line.rstrip().split("\t")
    b[0] = "chr1"
    sys.stdout.write("\t".join(b)+"\n")
