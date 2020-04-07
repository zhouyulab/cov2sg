#!/usr/bin/env python
"""
Get (DL, DR, AL, AR) for junctions site
      D       A
------         -----
   DL   DR  AL  AR
   +    +   -    -


"""
import os
import sys

if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.stderr.write("Usage: program f_jd armsize out.bed\n")
        sys.exit(1)

    f_jd, armsize, f_bed = sys.argv[1:4]
    armsize = int(armsize)
    assert os.path.exists(f_jd)
    assert armsize > 0
    foh = open(f_bed, "w")
    for line in open(f_jd):
        b = line.rstrip().split("\t")
        gid, jid, chrom, start, end, _, score, strand = b[:8]
        if strand != "+":
            continue

        s = int(start)
        e = int(end)
        arms = [
            ('DL', s-armsize, s, '+'),
            ('AL', e-armsize, e, '-'),
            ('DR', s, s+armsize, '+'),
            ('AR', e, e+armsize, '-'),
        ]
        for arm_name, arm_s, arm_e, arm_strand in arms:
            foh.write("\t".join(map(str, [
                chrom, arm_s, arm_e,
                "%d|%d|%s|%s|%s" % (s, e, arm_name, gid, jid), 0, arm_strand
            ]))+"\n")

    foh.close()
