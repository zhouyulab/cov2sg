import pysam
from plastid import GenomicSegment


class SpliceJuncDB:
    """Set of splicing junctions"""
    def __init__(self, sjdb_tbi):
        self.sjdb = pysam.TabixFile(sjdb_tbi)

    def fetch(self, chrom, start, end, strand=None):
        for hit in self.sjdb.fetch(chrom, start, end):
            seg = self.as_seg(hit)
            if strand is None:
                yield seg
            elif seg.strand == strand:
                yield seg

    def findseg(self, seg):
        assert isinstance(seg, GenomicSegment)
        for s in self.fetch(seg.chrom, seg.start, seg.end, seg.strand):
            yield s

    @staticmethod
    def as_seg(hit):
        chrom, start, end, strand = hit.split("\t")
        return GenomicSegment(chrom, int(start), int(end), strand)

    def has_splicejunc(self, seg, prange=1):
        for hit in self.findseg(seg):
            if (abs(hit.start - seg.start) <= prange and
                    abs(hit.end - seg.end) <= prange):
                return True
        return False
