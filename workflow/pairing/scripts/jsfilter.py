#!/usr/bin/env python
import os
import sys


if len(sys.argv) < 2:
    sys.stderr.write("Usage: program js1 js2 > js.csv\n")
    sys.exit(1)

fjs1, fjs2 = sys.argv[1:3]
assert os.path.exists(fjs1), fjs1
assert os.path.exists(fjs2), fjs2

sig2js = {}
for line in open(fjs1):
    b = line.rstrip().split("\t")
    gid, jid, chrom, start, end, _, _, strand = b[:8]
    sig = (int(start), int(end), strand)
    assert sig not in sig2js, sig
    sig2js[sig] = {
        'gid': gid,
        'jid': jid,
        "ntot": int(b[-1]),
        "chrom": chrom,
        'sig': sig
    }


class JG:
    def __init__(self, gid):
        self.gid = gid
        self.jsli = []

    def add(self, js):
        self.jsli.append(js)

    def ntot(self):
        return sum(js['ntot'] for js in self.jsli)

    def core(self):
        totli = [js['ntot'] for js in self.jsli]
        idx = totli.index(max(totli))
        return self.jsli[idx]


gid2jg = {}
for line in open(fjs2):
    if line.startswith("Start"):
        continue

    start, end, ngs1, ngs2, nano1, nano2 = list(map(int, line.rstrip().split("\t")))
    sig = (start, end, "+")
    assert sig in sig2js, sig
    js = sig2js[sig]
    gid = js['gid']
    if gid not in gid2jg:
        gid2jg[gid] = JG(gid)
    gid2jg[gid].add(js)


for jg in sorted(gid2jg.values(), key=lambda x: -x.ntot()):
    js = jg.core()
    name = "%s|%s" % (js["gid"], js["jid"])
    start, end, strand = js["sig"]
    sys.stdout.write("\t".join(map(str, [
        chrom, start, end, name, jg.ntot(), strand
    ]))+"\n")

