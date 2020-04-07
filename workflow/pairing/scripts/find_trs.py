#!/usr/bin/env python

import os
from Bio import SeqIO
from Bio.Seq import Seq

R = os.path.join(os.environ["HOME"], "gitee/zhouyulab/sgcov")
resdir = os.path.join(R, "results")
g = [r for r in SeqIO.parse(os.path.join(R, "genome/WIV04.fasta"), "fasta")][0]
print(len(g))
gs = g.seq

trs = Seq("TCTAAACGAACTTTAA")
idx = gs.find(trs)
while idx >= 0:
    print(idx)
    idx = gs.find(trs, idx+1)

trsdir = os.path.join(resdir, "trs")
if not os.path.exists(trsdir):
    os.mkdir(trsdir)

with open(os.path.join(trsdir, "trs.fa"), "w") as foh:
    foh.write(">trs\n%s\n" % trs)

