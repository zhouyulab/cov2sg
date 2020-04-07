#!/usr/bin/env python
"""
BAM file readers and grouping iterators
"""
from collections import defaultdict

def ibam_chrom(bamin):
    """Iterator of reads from sorted BAM by chrom"""
    for chrom in bamin.references:
        yield chrom, bamin.fetch(chrom)


def ibam_chromStart(bamin):
    """
    Iterator of reads grouped by signature(chrom, start)
    yield (sig, [reads])
    bamin is readed from stdin from bam files sorted by samtools sort
    """
    for chrom, ir in ibam_chrom(bamin):
        psig = None # position signature
        reads = []
        for read in ir:
            if read.is_unmapped:
                continue

            sig = (chrom, read.pos)
            if psig is None:
                psig = sig
            elif psig != sig:
                yield (psig, reads)
                reads = [] # reset
                psig = sig

            reads.append(read)

        if reads:
            yield (psig, reads)


def igroup_strand(reads):
    """
    Group given reads by strand
    yield (strand, [reads])
    """
    fwdreads = []
    revreads = []
    for read in reads:
        if read.is_reverse:
            revreads.append(read)
        else:
            fwdreads.append(read)
    if fwdreads:
        yield ('+', fwdreads)

    if revreads:
        yield ('-', revreads)


def igroup_end(reads):
    """
    Group given reads by end position
    yield (end, [reads])
    """
    end2reads = defaultdict(list)
    for read in reads:
        end2reads[read.aend].append(read)

    for end in sorted(end2reads):
        yield (end, end2reads[end])


def ibam_chromStartEndStrand(bamin):
    """
    Iterator of reads grouped by signature(chrom, start, end, strand)
    yield (sig, [reads])
    bamin is readed from stdin from bam files sorted by samtools sort
    """
    for (chrom, start), reads1 in ibam_chromStart(bamin):
        for end, reads2 in igroup_end(reads1):
            for strand, reads3 in igroup_strand(reads2):
                yield ((chrom, start, end, strand), reads3)


def reads2qnames(reads):
    """Return qname list from given reads"""
    return [read.qname for read in reads]


def igroup_barcode(reads, bclen=4):
    """([], int) -> yield (str, [])
    Iterator of reads grouped by signature(barcode)
    barcode is assumed to be at the end of the read qname
    """
    bc2reads = defaultdict(list)
    for read in reads:
        name = read.qname
        bc = name[(len(name)-bclen):]
        bc2reads[bc].append(read)

    for bc, reads in bc2reads.items():
        yield (bc, reads)


def ibampe(bamin):
    """Iterator of proper-pair PE-reads from sorted BAM file"""
    for chrom, ir in ibam_chrom(bamin):
        rnli = []
        rn2pe = defaultdict(dict)
        for r in ir:
            if r.is_paired and not r.is_secondary:
                if r.qname not in rn2pe:
                    rnli.append(r.qname)
                if r.is_read1:
                    rn2pe[r.qname]['R1'] = r
                elif r.is_read2:
                    rn2pe[r.qname]['R2'] = r
        yield chrom, [rn2pe[rn] for rn in rnli if rn in rn2pe]


def ibampe_chromStartEndR1Strand(bamin):
    """Iterator of PE-reads by (chrom, start, end, R1Strand)"""
    for chrom, rpeli in ibampe(bamin):
        sig2rli = defaultdict(list)
        for pe in rpeli:
            if not ('R1' in pe and 'R2' in pe):
                continue
            start = min(pe['R1'].reference_start, pe['R2'].reference_start)
            end = max(pe['R1'].reference_end, pe['R2'].reference_end)
            strand = '+'
            if pe['R1'].is_reverse:
                strand = '-'

            sig = (start, end, strand)
            sig2rli[sig].append(pe)

        for sig in sorted(sig2rli):
            yield (chrom, sig[0], sig[1], sig[2]), sig2rli[sig]

