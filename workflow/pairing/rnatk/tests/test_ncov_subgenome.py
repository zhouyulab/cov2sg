import pysam
import pytest
from plastid import GenomicSegment
from rnatk.ncov.subgenome import Jnc, JncAssembler, JncGroup, JncSite, SegRead


@pytest.fixture()
def bam():
    SegRead.LIBTYPE = None
    return pysam.AlignmentFile("tests/data/trsread.bam", "rb")


@pytest.fixture(scope="module")
def rex():
    SegRead.LIBTYPE = None
    bam = pysam.AlignmentFile("tests/data/trsread.bam", "rb")
    return [r for r in bam.fetch() if r.qname == "1"][0]


class TestSegRead:
    def test_name(self, rex):
        sr = SegRead(rex)
        assert sr.name == '1'

    def test_sgc(self, rex):
        sr = SegRead(rex)
        assert str(sr.sgc) == "MN996528:105-130^330-355(+)"

    def test_g_junc(self, rex):
        sr = SegRead(rex)
        juncs = list(sr.g_jnc())
        assert len(juncs) == 1

    def test_ref_loc2query_loc(self, rex):
        sr = SegRead(rex)
        assert sr.ref_loc2query_loc(105, 110) == (0, 5)
        assert sr.ref_loc2query_loc(330, 335) == (25, 30)
        assert sr.ref_loc2query_loc(125, 335) == (20, 30)

    def test_fetch_sequence(self, rex):
        sr = SegRead(rex)
        s = "CATGCTTAGTGCACTCACGCAGTATGGTTCGCGACGTGCTCGTACGTGGC"
        assert sr.fetch_sequence(105, 110) == s[0:5]
        assert sr.fetch_sequence(330, 335) == s[25:30]
        assert sr.fetch_sequence(125, 335) == s[20:30]


class TestJnc:
    def test_lft(self, rex):
        sr = SegRead(rex)
        juncs = list(sr.g_jnc())
        assert str(juncs[0].lft) == "MN996528:105-130(+)"

    def test_rgt(self, rex):
        sr = SegRead(rex)
        juncs = list(sr.g_jnc())
        assert str(juncs[0].rgt) == "MN996528:330-355(+)"

    def test_has_minlen(self, rex):
        sr = SegRead(rex)
        juncs = list(sr.g_jnc())
        assert juncs[0].has_minlen() is True
        tmp = Jnc.ARM_MINLEN
        Jnc.ARM_MINLEN = 30
        assert juncs[0].has_minlen() is False
        Jnc.ARM_MINLEN = tmp

    def test_length(self, rex):
        sr = SegRead(rex)
        juncs = list(sr.g_jnc())
        assert juncs[0].length == 200


class TestJncSite:
    def test_jncsite1(self, bam):
        ja = JncAssembler(bam, window=1, minevi=1)
        ja.get_jncs("MN996528", 1, 500, "+")
        assert len(ja.jsdi) == 2
        gs1 = GenomicSegment("MN996528", 130, 330, "+")
        js1 = ja.jsdi[gs1]
        assert js1.ntot == 4
        assert js1.nevi == 3
        gs2 = GenomicSegment("MN996528", 130, 332, "+")
        js2 = ja.jsdi[gs2]
        assert js2.ntot == 1
        assert js2.nevi == 1
        assert js2.as_bed() == ['MN996528', '130', '332', 'MN996528:130-131^331-332(+)', '0', '+',
                                '130', '332', '0,0,0', '2', '1,1,', '0,201,', '1', '1']


class TestJncGroup:
    def test_jncgroup1(self, bam):
        ja = JncAssembler(bam, window=1, minevi=1)
        ja.get_jncs("MN996528", 1, 500, "+")
        ja.cluster()
        assert len(ja.jgli) == 1
        jg = ja.jgli[0]
        assert jg.get_ntot() == 5
        assert jg.get_nevi() == 4
        assert str(jg.corejs) == "MN996528\t130\t330\t+\t4\t3"
        jg.merge()
        jg.set_name("JncGroup")
        assert jg.as_bed() == list(map(str, [
            "MN996528", 130, 332, "JncGroup", 0, "+", 130, 130, "255,0,0", 3, "1,1,1,", "0,199,201,", 4, 5]))
        jsli = [tuple(js[:5] + js[-2:]) for js in jg.jncsites()]
        assert jsli == [('JncGroup', '1', 'MN996528', '130', '330', '3', '4'),
                        ('JncGroup', '2', 'MN996528', '130', '332', '1', '1')]


class TestJncAssembler:
    def test_get_jncs(self, bam):
        ja = JncAssembler(bam, window=1, minevi=1)
        assert ja.numjs() == 0
        ja.get_jncs("MN996528", 1, 500, "+")
        assert ja.numjs() == 2

    def test_sorted(self, bam):
        ja = JncAssembler(bam, window=0, minevi=1)
        ja.get_jncs("MN996528", 1, 500, "+")
        jsli = [js for _, js in ja.sorted()]
        assert len(jsli) >= 2
        assert (jsli[0].nevi, jsli[0].ntot) >= (jsli[1].nevi, jsli[1].ntot)

    def test_cluster1(self, bam):
        ja = JncAssembler(bam, window=0, minevi=1)
        ja.get_jncs("MN996528", 1, 500, "+")
        ja.cluster()
        assert ja.numjg() == 2

    def test_cluster2(self, bam):
        ja = JncAssembler(bam, window=1, minevi=1)
        ja.get_jncs("MN996528", 1, 500, "+")
        ja.cluster()
        assert ja.numjg() == 1

    def test_merge1(self, bam):
        ja = JncAssembler(bam, window=1, minevi=1)
        ja.get_jncs("MN996528", 1, 500, "+")
        ja.cluster()
        ja.merge()
        assert str(ja.jgli[0].sc) == "MN996528:130-131^329-330^331-332(+)"

    def test_assemble(self, bam):
        ja = JncAssembler(bam, window=5, minevi=1)
        ja.assemble("MN996528", 1, 2000, "+")
        assert ja.numjg() == 2
        assert str(ja.jgli[0]) == "\t".join(map(str, ["MN996528", 130, 330, "+", 4, 3, 5, 4]))
        assert str(ja.jgli[1]) == "\t".join(map(str, ["MN996528", 130, 1130, "+", 1, 1, 2, 2]))

    def test_i_chromsize(self, bam):
        ja = JncAssembler(bam, window=5, minevi=1)
        csli = list(ja.i_chromsize())
        assert csli == [('MN996528', 29891)]

    def test_assemble_checkseq1(self, bam):
        genome = "tests/data/WIV04.fasta"
        ja = JncAssembler(bam, window=5, minevi=1, genome=None)
        ja.assemble("MN996528", 800, 4000, "+")
        assert ja.numjg() == 1
        assert ja.jgli[0].get_ntot() == 3

    def test_assemble_checkseq2(self, bam):
        genome = "tests/data/WIV04.fasta"
        Jnc.ARM_MINLEN = 20
        ja = JncAssembler(bam, window=5, minevi=1, genome=genome)
        ja.assemble("MN996528", 800, 4000, "+")
        assert ja.numjg() == 1
        assert ja.jgli[0].get_ntot() == 1
        Jnc.ARM_MINLEN = 8

    def test_assemble_checkseq3(self, bam):
        genome = "tests/data/WIV04.fasta"
        Jnc.ARM_MINLEN = 11
        ja = JncAssembler(bam, window=5, minevi=1, genome=genome)
        ja.assemble("MN996528", 800, 4000, "+")
        assert ja.numjg() == 1
        assert ja.jgli[0].get_ntot() == 2
        Jnc.ARM_MINLEN = 8

    def test_assemble_checkseq4(self, bam):
        genome = "tests/data/WIV04.fasta"
        Jnc.ARM_MINLEN = 10
        ja = JncAssembler(bam, window=5, minevi=1, genome=genome)
        ja.assemble("MN996528", 800, 4000, "+")
        assert ja.numjg() == 1
        assert ja.jgli[0].get_ntot() == 3
        Jnc.ARM_MINLEN = 8
