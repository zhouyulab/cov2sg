"""
Fastq format
"""


class Fastq(object):
    """
    Fastq record
    """
    def __init__(self, name, seq, qual):
        self.name = name
        self.seq = seq
        self.qual = qual
    
    def __str__(self):
        return '\n'.join(['@%s' % str(self.name), self.seq, 
            '+', self.qual])
    
    def as_simple(self):
        return '\n'.join(['@%s' % str(self.name), self.seq, '+', self.qual])
    
    def get_name(self):
        return self.name
    
    def clean_name(self):
        self.name = self.name.split(" ")[0]
    
    def get_seq(self):
        return self.seq
    
    def get_qual(self):
        return self.qual
    
    def get_ascii(self):
        return map(ord, self.qual)
    
    def get_length(self):
        """Return length"""
        return len(self.seq)

    def set_name(self, name):
        """Reset name"""
        self.name = name

    
def reader_fastq(infile):
    """Generator of Fastq object from given file"""
    i = 0
    name = None
    seq = None
    qual = None
    for line in open(infile):
        i += 1
        curr_line = line.strip()
        if i % 4 == 1:
            name = curr_line[1:]
        elif i % 4 == 2:
            seq = curr_line
        elif i % 4 == 0:
            qual = curr_line
            yield Fastq(name, seq, qual)
            

def readfq(fp): # this is a generator function
    """https://github.com/lh3/readfq/blob/master/readfq.py"""
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break
