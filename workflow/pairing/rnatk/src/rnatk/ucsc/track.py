"""
UCSC track utilites
"""
import sys
import subprocess
from bx.tabular.io import Comment
from rnatk.ios.bed import BedReader, Bed


class TrackMaker(object):
    """
    Making track
    """
    def __init__(self, rgb_fwd='255,0,0', rgb_rev='0,0,255'):
        self.strand2rgb = {'+':rgb_fwd, '-':rgb_rev, }

    def get_thick(self, b):
        """Get thick start and end position"""
        return b.get_start(), b.get_end()

    def gr_bed12(self, infile, keep_comment=False):
        """Generator of BED12 lines"""
        blockCount = 1
        blockStarts = '0,'
        for b in BedReader(open(infile), return_comments=keep_comment,
            comment_lines_startswith = ["track "]):
            if isinstance(b, Comment):
                yield b.line
            if isinstance(b, Bed):
                thickStart, thickEnd = self.get_thick(b)
                blockSizes = '%d,' % (b.get_end() - b.get_start(), )
                try:
                    score = int(float(b.get_score()))
                except(ValueError):
                    score = 0

                yield '\t'.join(map(str,
                    [b.get_chrom(), b.get_start(), b.get_end(),
                    b.get_name(), score, b.get_strand(),
                    thickStart, thickEnd, self.strand2rgb[b.get_strand()],
                    blockCount, blockSizes, blockStarts]))

    def bedTo12(self, infile, outfile, keep_comment=False, sort=False):
        """
        To bed format of 12 columns
        """
        foh = open(outfile, "w")
        for b in self.gr_bed12(infile, keep_comment):
            foh.write(b+"\n")
        foh.close()

        if sort:
            retcode = 0
            try:
                retcode = subprocess.call(['bedSort', outfile, outfile])
            except OSError:
                sys.stderr.write('No bedSort program from kent source\n')
                return False

            if retcode != 0:
                sys.stderr.write('bedSort not executed successfully\n')
                return False

        return True

    def bed6To12(self, infile, outfile, keep_comment=False, sort=False):
        """Keep it"""
        self.bedTo12(infile, outfile, keep_comment, sort)


class TrackBed8(TrackMaker):
    """Bed8 format"""
    def get_thick(self, b):
        """Get thick start and end position"""
        return b.get_thickStart(), b.get_thickEnd()


