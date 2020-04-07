#!/usr/bin/env python
"""Module for RNAhybrid
"""
from Bio import SeqIO
import os
import sys
import argparse
import logging
_logger = logging.getLogger(__name__)


class Target:
    """
    Target from using compact mode of RNAhybrid
    """
    def __init__(self, compact_target):
        fields = compact_target.rstrip("\n").split(":")
        assert len(fields) == 11
        lens = [len(fields[i]) for i in range(7, 11)]
        assert len(set(lens)) == 1, 'not correct pairing lines: %s' % "".join(map(str, lens))
        self.fields = fields

    def __str__(self):
        return "\t".join([self.query + self.target])

    @property
    def tname(self):
        return self.fields[0]

    @property
    def qname(self):
        return self.fields[2]

    @property
    def mfe(self):
        return float(self.fields[4])

    @property
    def pvalue(self):
        return float(self.fields[5])

    @property
    def position(self):
        return int(self.fields[6])

    @property
    def pairing(self):
        return self.fields[7:11]

    @property
    def pairing_rc(self):
        return ["".join(reversed(line)) for line in reversed(self.pairing)]

    @property
    def target(self):
        """Return target 5'->3' sequence and base pairing"""
        t = ""
        b = ""
        unpaired = self.fields[7]
        paired = self.fields[8]
        assert len(unpaired) == len(paired)
        i = 0
        while i < len(paired):
            if unpaired[i] != " " or paired[i] != " ":
                if unpaired[i] != " ":
                    t += unpaired[i]
                    b += "."
                else:
                    t += paired[i]
                    b += "*"
            i += 1

        return t, b

    @property
    def query(self):
        """Return query 5'->3' sequence and base pairing"""
        q = ""
        b = ""
        paired = self.fields[9]
        unpaired = self.fields[10]
        assert len(unpaired) == len(paired)
        i = 0
        while i < len(paired):
            if paired[i] != " " or unpaired[i] != " ":
                if paired[i] != " ":
                    q += paired[i]
                    b += "*"
                else:
                    q += unpaired[i]
                    b += "."
            i += 1

        # query in the bottom as 3'->5'
        return "".join(reversed(q)), "".join(reversed(b))

    def target_full(self, tgt_full):
        tgt_full = tgt_full.upper().replace("T", "U")
        t, s = self.target
        assert tgt_full[(self.position-1):(self.position-1+len(t))] == t, "full target sequence not contain hit target"
        s_full = ["."] * len(tgt_full)
        for i in range(len(s)):
            s_full[i+self.position-1] = s[i]

        return tgt_full, "".join(s_full)

    def as_cofold(self, tgt_full):
        qs, qb = self.query
        ts, tb = self.target_full(tgt_full)
        cos = qs + "&" + ts
        cob = qb.replace("*", "(") + "&" + tb.replace("*", ")")
        return cos, cob


def parse_args(args):
    """Parse command line parameters"""
    parser = argparse.ArgumentParser(
        description="Parse parameters for %s" % __file__)
    parser.add_argument(
        "hyb_in",
        type=str,
        help="Input RNAhybrid -c output file",
        metavar="hybrid.csv")
    parser.add_argument(
        "pair_fa",
        type=str,
        help="Fasta file for paired query/target sequences",
        metavar="pair.fa")
    parser.add_argument(
        "cofold_out",
        type=str,
        help="Output sequence and structure in RNAcofold format",
        metavar="cofold.csv")
    parser.add_argument(
        '-v',
        '--verbose',
        dest="loglevel",
        help="set loglevel to INFO",
        action='store_const',
        const=logging.INFO)
    parser.add_argument(
        '-vv',
        '--very-verbose',
        dest="loglevel",
        help="set loglevel to DEBUG",
        action='store_const',
        const=logging.DEBUG)
    return parser.parse_args(args)


def setup_logging(loglevel):
    """Setup basic logging
    Args:
      loglevel (int): minimum loglevel for emitting messages
    """
    logformat = "[%(asctime)s] %(levelname)s:%(name)s:%(message)s"
    logging.basicConfig(level=loglevel, stream=sys.stdout,
                        format=logformat, datefmt="%Y-%m-%d %H:%M:%S")


def main(args):
    """Main entry point allowing external calls"""
    args = parse_args(args)
    setup_logging(args.loglevel)
    _logger.debug("Starting reading ...")
    assert os.path.exists(args.pair_fa)
    sid2qt = {}
    for r in SeqIO.parse(args.pair_fa, "fasta"):
        query, target = str(r.seq).split("&")
        sid2qt[r.id] = {'query': query, "target": target}

    with open(args.cofold_out, "w") as foh:
        i = 0
        for line in open(args.hyb_in):
            i += 1
            t = Target(line)
            gid, jid, _ = t.tname.split("_")
            sig = "%s|%s" % (gid, jid)
            q_full = sid2qt[sig]["query"].replace("T", "U")
            t_full = sid2qt[sig]["target"]
            q, _ = t.query
            assert q == q_full, "%s != %s" % (q, q_full)
            jseq, jstruct = t.as_cofold(t_full)
            foh.write(" ".join(map(str, [
                i, sig, jseq, jstruct, t.mfe
            ]))+"\n")

    _logger.info("Ends.")


def run():
    """Entry point for console_scripts"""
    main(sys.argv[1:])


if __name__ == '__main__':
    run()
