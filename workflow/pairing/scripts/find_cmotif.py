#!/usr/bin/env python
"""
Find core motifs (ACGAAC/AAGAAC) JV2004
"""
import os
import sys
from Bio import SeqIO


def find_all(g, motif):
    idx = 0
    while True:
        idx = g.find(motif, idx)
        if idx == -1:
            return
        yield idx
        idx += 1


if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.stderr.write("Usage: program in.fa out.bed\n")
        sys.exit(1)

    fa_in, f_bed = sys.argv[1:3]
    assert os.path.exists(fa_in)
    r = [r for r in SeqIO.parse(fa_in, "fasta")][0]
    gfwd = r.seq
    grev = r.seq.reverse_complement()
    gsize = len(r.seq)
    motifs = ('ACGAAC', "AAGAAC")
    with open(f_bed, "w") as foh:
        for m in motifs:
            for strand in ('+', '-'):
                if strand == '+':
                    gs = gfwd
                else:
                    gs = grev
                for start in find_all(gs, m):
                    if strand == '-':
                        start = gsize - len(m) - start
                    foh.write("\t".join(map(str, [
                        r.id, start, start+len(m), m, 0, strand]))+"\n")
