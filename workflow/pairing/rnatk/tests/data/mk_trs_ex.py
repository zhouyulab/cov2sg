#!/usr/bin/env python

import os
from Bio import SeqIO

fagenome = "WIV04.fasta"
g = [r for r in SeqIO.parse(fagenome, "fasta")][0]
gs = g.seq

# Construct lead sequence and body sequence test dataset
"""
    ---         --- 
  -----         -----
  -----         -----
 ------         ------
-------         ----- 
  -----                      ------
  -----                      ------
"""
RL = 30
GAP = 200
reads = [
    (105, RL-5, GAP, RL-5),
    (100, RL, GAP, RL),
    (100, RL, GAP, RL),
    (99,  RL+1, GAP, RL+1),
    (98,  RL+2, GAP+2, RL),
    (100, RL, GAP*5, RL),
    (100, RL, GAP*5+5, RL),
    (2100, RL, GAP * 5, RL, None, None),
    (2100, RL, GAP * 5, RL, -20, None),
    (2100, RL, GAP * 5, RL, None, 10),
]

with open("trsread.fa", "w") as foh:
    i = 0
    for r in reads:
        if len(r) == 4:
            start, lsgl, gapl, rsgl = r
        elif len(r) == 6:
            start, lsgl, gapl, rsgl, mlft, mrgt = r
        i += 1
        lsg = gs[start:start+lsgl]
        rstart = start + lsgl + gapl
        rsg = gs[rstart:rstart+rsgl]

        if len(r) == 6 and mlft is not None:
            lsgli = list(lsg)
            lsgli[mlft] = lsg.complement()[mlft]
            lsg = "".join(lsgli)

        if len(r) == 6 and mrgt is not None:
            rsgli = list(rsg)
            rsgli[mrgt] = rsg.complement()[mrgt]
            rsg = "".join(rsgli)

        jcr = lsg + rsg
        foh.write(">%d\n%s\n" % (i, jcr))
