#!/usr/bin/env python
import os
import sys


if len(sys.argv) < 3:
    sys.stderr.write("Usage: program jsseq.csv outfile\n")
    sys.exit(1)

fcsv, outfile = sys.argv[1:3]
assert os.path.exists(fcsv), fcsv

outdir = outfile + ".rnahyb"
if not os.path.exists(outdir):
    os.mkdir(outdir)

foh_csv = open(outfile, "w")
for line in open(fcsv):
    b = line.rstrip().split("\t")
    gid, jid = b[3].split("|")
    qid = "%s_%s_q" % (gid, jid)
    tid = "%s_%s_t" % (gid, jid)
    qfa = os.path.join(outdir, "%s.fa" % qid)
    tfa = os.path.join(outdir, "%s.fa" % tid)

    foh_q = open(qfa, "w")
    foh_q.write(">%s\n%s\n" % (qid, b[8]))
    foh_q.close()

    foh_t = open(tfa, "w")
    foh_t.write(">%s\n%s\n" % (tid, b[10]))
    foh_t.close()

    fhyb = os.path.join(outdir, "%s_%s.csv" % (gid, jid))
    cmdline = "RNAhybrid -b 1 -c -s 3utr_human -t %s -q %s > %s" % (tfa, qfa, fhyb)
    os.system(cmdline)
    os.unlink(qfa)
    os.unlink(tfa)

    foh_csv.write(open(fhyb).readline())

foh_csv.close()
