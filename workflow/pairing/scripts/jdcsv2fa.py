#!/usr/bin/env python
import os
import sys


if len(sys.argv) < 2:
    sys.stderr.write("Usage: program jd.csv > js.csv\n")
    sys.exit(1)

fcsv = sys.argv[1]
assert os.path.exists(fcsv), fcsv
for line in open(fcsv):
    b = line.rstrip().split("\t")
    rid = b[3]
    seq = b[8] + "&" + b[10]
    sys.stdout.write(">%s\n%s\n" % (rid, seq))
