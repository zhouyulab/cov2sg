from bx.tabular.io import ParseError, TableRow, TableReader, Header
import sys


def read_bed(fname):
    """Simple function for BED6"""
    fhd = open(fname)
    for line in fhd:
        if line.startswith('#'):
            continue
        fields = line.rstrip( "\r\n" ).split( "\t" )
        yield (fields[0], int(fields[1]), int(fields[2]), 
               fields[3], fields[4], fields[5])
    fhd.close()


class BedFieldNumberError(ParseError):
    """Exception"""
    pass


class BedFieldFormatError(ParseError):
    """Exception"""
    pass


class Bed(TableRow):
    """
    A Bed record
    """
    def __init__(self, fields, reader=None):
        if reader != None:
            TableRow.__init__(self, reader, fields)
        else:
            self.fields = fields

        self.nfields = len(fields)
        if self.nfields < 3:
            raise BedFieldNumberError("Not enough fields, at least 3")

    def __str__( self ):
        return "\t".join(map(str, self.fields))
    
    def __cmp__(self, other):
        return cmp(self.fields, other.fields)
    
    def get_nfields(self):
        """number of fields"""
        return self.nfields
    
    def get_start(self):
        """start pos"""
        return int(float(self.__getitem__('chromStart')))
    
    def get_end(self):
        """end pos"""
        return int(float(self.__getitem__('chromEnd')))
    
    def get_chrom(self):
        """get chrom"""
        return self.__getitem__('chrom')
    
    def get_strand(self):
        """get strand"""
        return self.__getitem__('strand')
    
    def get_name(self):
        return self.__getitem__('name')
    
    def get_score(self):
        return self.__getitem__('score')
    
    def get_thickStart(self):
        return self.__getitem__('thickStart')
    
    def get_thickEnd(self):
        return self.__getitem__('thickEnd')
    
    def get_blockCount(self):
        return self.__getitem__('blockCount')
    
    def get_blockSizes(self):
        return self.__getitem__('blockSizes')
    
    def get_blockStarts(self):
        return self.__getitem__('blockStarts')
    
    def copy(self):
        """
        Create new Bed record by copying self
        """
        return Bed(list(self.fields), self.reader)
    
    def get_length(self):
        """Return length"""
        return self.get_end() - self.get_start()
    
    def find_middle(self):
        """
        Find middle point position
        """
        return self.get_start() + self.get_length()/2
    
    def find_midpoint_interval(self):
        """
        Return 1-length middle point interval
        """
        mid = self.find_middle()
        if self.get_strand() == '+':
            s = mid
            e = s + 1
        else:
            e = mid
            s = mid - 1
        return (s, e)
        
    def extend_twosides(self, extlen=100, by=True, chrom2size=None):
        """Extend to two sides"""
        return self.extend('twosides', extlen, by, chrom2size)

    def extend_3end(self, extlen=100, by=True, chrom2size=None):
        """Extend to 3'end"""
        if self.get_strand() == '+':
            return self.extend('right', extlen, by, chrom2size)
        else:
            return self.extend('left', extlen, by, chrom2size)

    def extend_5end(self, extlen=100, by=True, chrom2size=None):
        """Extend to 3'end"""
        if self.get_strand() == '+':
            return self.extend('left', extlen, by, chrom2size)
        else:
            return self.extend('right', extlen, by, chrom2size)    
    
    def extend(self, direction, extlen=100, by=True, chrom2size=None):
        """
        Return new Bed record with extension
        direction: twosides|left|right
        extlen: extended length
        by: True, add extlen according to direction; 
            False, only extend when current length < extlen;
        chrom2size: if True, chop the extended record 
        """
        assert direction in ['twosides', 'left', 'right']
        bed_ext = self.copy()
        if by:
            if direction == 'twosides':
                bed_ext.fields[1] = str(self.get_start() - extlen)
                bed_ext.fields[2] = str(self.get_end() + extlen)
            elif direction == 'left':
                bed_ext.fields[1] = str(self.get_start() - extlen)
            else:
                bed_ext.fields[2] = str(self.get_end() + extlen)
        else:
            if self.get_length() >= extlen:
                return bed_ext
            else:
                if direction == 'twosides':
                    middle = self.find_middle()
                    leftlen = int(extlen / 2)
                    bed_ext.fields[1] = str(middle - leftlen)
                    bed_ext.fields[2] = str(middle + extlen - leftlen)
                elif direction == 'left':
                    bed_ext.fields[1] = str(self.get_end() - extlen)
                else:
                    bed_ext.fields[2] = str(self.get_start() + extlen)
                    
        if chrom2size:
            bed_ext.chop(chrom2size)
        return bed_ext
    
    def chop(self, chrom2size):
        """Chop start, end so as to not pass chromStart, chromEnd"""
        chrom = self.get_chrom()
        assert chrom in chrom2size        
        self.fields[1] = str(max(0, self.get_start()))
        self.fields[2] = str(min(chrom2size[chrom], self.get_end()))



class BedReader(TableReader):
    """
    Reader for Bed format
    """
    fieldNames = ['chrom', 'chromStart', 'chromEnd', 'name', 'score',
        'strand', 'thickStart', 'thickEnd', 'itemRgb',
        'blockCount', 'blockSizes', 'blockStarts']
    
    def __init__(self, fhd, discard_first_column=False, return_header=False,
                 return_comments=False, force_header=None,
                 comment_lines_startswith = ["#", "track ", "browser"]):
        
        TableReader.__init__(self, fhd, return_header, return_comments, 
                             force_header, comment_lines_startswith)
        self.discard_first_column = discard_first_column
        if not self.header:
            self.header = Header(BedReader.fieldNames)

    def parse_row(self, line):
        """return Bed object"""
        fields = line.split("\t")
        if self.discard_first_column:
            fields.pop(0)
        return Bed(fields, self) #self as argument reader

    
def split_by_strand(infile, outfile_f, outfile_r, clean=None):
    """
    Split Bed files by strand.
    Generate two files: outfile_f, outfile_r
    clean: whether remove other fields except chrom, chromStart, chromEnd.
    """
    outhandle_f = open(outfile_f, "w")
    outhandle_r = open(outfile_r, "w")
    
    breader = BedReader(open(infile))
    for b in breader:
        assert isinstance(b, Bed)
        if b.get_strand() == '+':
            if not clean:
                outhandle_f.write(str(b)+'\n')
            else:
                outhandle_f.write("\t".join(
                    [b.get_chrom(), str(b.get_start()), str(b.get_end())])+'\n')
        elif b.get_strand() == '-':
            if not clean:
                outhandle_r.write(str(b)+'\n')
            else:
                outhandle_r.write("\t".join(
                    [b.get_chrom(), str(b.get_start()), str(b.get_end())])+'\n')
    
    outhandle_f.close()
    outhandle_r.close()
    
    sys.stderr.write("Two files generated: %s, %s\n" % (outfile_f, outfile_r))
