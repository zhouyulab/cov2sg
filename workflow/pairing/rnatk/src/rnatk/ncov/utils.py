
class Cofold:
    def __init__(self, b):
        seq_num, seq_id, seq, mfe_struct, mfe, p_struct, p_mfe, p_freq, deltag = list(
            map(str.strip, b))
        self.id = seq_id
        self.seq = seq.strip('"')
        self.mfe_struct = mfe_struct.strip('"')
        self.mfe = float(mfe)
        self.p_struct = p_struct.strip('"')
        self.p_mfe = float(p_mfe)
        self.p_freq = float(p_freq)
        self.deltag = float(deltag)

    def __str__(self):
        return "\n".join(self.info)

    @property
    def info(self):
        return [
            "\t".join(map(str, [self.id, self.seq, self.seq,  self.mfe])),
            "\t".join(map(str, [self.id, self.mfe_struct, self.p_struct, self.p_mfe]))]


def read_cofold(fcofoldp):
    d = {}
    i = 0
    for line in open(fcofoldp, encoding="ISO-8859-1"):
        i += 1
        if i == 1:
            idx = line.index("mfe ")
            line = line[idx+4:]
        b = line.rstrip().split(" ")
        b = [v for v in b if v != ""]
        cf = Cofold(b)
        d[cf.id] = cf
    return d


class JS:
    def __init__(self, b):
        chrom, start, end, name, _, strand, numevi, ntot, fs1, fs2, fs2rc = b
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.id = name
        self.strand = strand
        self.numevi = numevi
        self.ntot = ntot
        self.fs1 = fs1
        self.fs2 = fs2
        self.fs2rc = fs2rc

    @property
    def sig(self):
        return self.start, self.end

    def __str__(self):
        return "\t".join(map(str, [self.chrom, self.start, self.end, self.id, self.numevi, self.ntot]))


def read_js(fjs):
    assert os.path.exists(fjs), fjs
    d = {}
    for line in open(fjs):
        b = line.rstrip().split("\t")
        js = JS(b)
        d[js.id] = js
    return d


class JncTwoSide:
    def __init__(self, jslft, cflft, jsrgt, cfrgt):
        assert jslft.sig == jsrgt.sig
        assert jslft.id == jsrgt.id
        self.jslft = jslft
        self.cflft = cflft
        self.jsrgt = jsrgt
        self.cfrgt = cfrgt

    @property
    def sig(self):
        return self.jslft.sig

    def format_struct(self, side):
        if side == 'lft':
            return ["LFT\t" + v for v in self.cflft.info]
        elif side == 'rgt':
            return ["RGT\t" + v for v in self.cfrgt.info]
        else:
            return None

    @property
    def info(self):
        return "\n".join(map(str, [
            self.jslft] + self.format_struct('lft') + self.format_struct('rgt')))+"\n"
