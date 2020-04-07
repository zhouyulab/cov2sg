#!/usr/bin/env python
"""
Classes related to Interval, Feature
"""


class Interval:
    """
    Basic interval, any object with start and end properties will work as well
    Originally from bx-python intervals.intersection.
    """
    def __init__(self, start, end, value=None):
        assert start <= end, "start must be less than end"
        self.start = start
        self.end = end
        self.value = value

    def __eq__(self, other):
        return self.start == other.start and self.end == other.end

    def __hash__(self):
        return hash((self.start, self.end))
    
    def __repr__(self):
        if self.value:
            return "IntervalWithValue(%d, %d, %s)" % \
                (self.start, self.end, repr(self.value))
        else:
            return "Interval(%d, %d)" % (self.start, self.end)

    def __str__(self):
        return "(%d, %d)" % (self.start, self.end)

    def get_start(self):
        """start pos"""
        return self.start
    
    def get_end(self):
        """end pos"""
        return self.end

    def get_len(self):
        """length of the interval"""
        return self.end - self.start

    def get_value(self):
        """associated value"""
        return self.value

    def set_value(self, value):
        """associated value"""
        self.value = value

    def overlaps(self, other):
        """whether overlaps with other interval"""
        return min(self.end, other.end) - max(self.start, other.start) > 0
    
    def write_pos(self):
        """return position string"""
        return '%d-%d' % (self.start, self.end)
    
    def get_midpoint(self):
        """get center coordinate"""
        return self.start + self.get_len() / 2
    
    def contains_pos(self, pos):
        return self.start <= pos < self.end
    

class Feature(Interval):
    """
    Genomic feature with chrom, strand
    """
    def __init__(self, chrom, start, end, strand='+', value=0):
        Interval.__init__(self, start, end, value)
        self.chrom = chrom
        self.strand = strand

    def __eq__(self, other):
        return (Interval.__eq__(self, other)
                and self.chrom == other.chrom
                and self.strand == other.strand)

    def __hash__(self):
        return hash((self.chrom, self.strand, self.start, self.end))

    def __str__(self):
        s = self.get_name() 
        if self.strand:
            s += ":" + self.strand
        if self.value:
            s += ":" + str(self.value)
        return s
        
    def get_chrom(self):
        """get chrom"""
        return self.chrom

    def get_name(self):
        """get name, default chrom:start-end"""
        return "%s:%d-%d" % (self.chrom, self.start, self.end)
    
    def get_bed(self):
        """get Bed format (6 cols) """
        return '\t'.join(map(str, [
            self.chrom, self.start, self.end,
            self.get_name(), self.value, self.strand]))
    
    def strand2symbol(self):
        """convert +/- to f/r"""
        if self.strand == '+':
            return 'f'
        elif self.strand == '-':
            return 'r'
        else:
            return ''

    def overlaps(self, other):
        """check whether overlap with other feature"""
        assert isinstance(other, Feature)
        if not (self.chrom == other.chrom and self.strand == other.strand):
            return False
        return Interval.overlaps(self, other)
