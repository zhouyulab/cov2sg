import pysam
import pytest
from rnatk.rna.hybrid import Target


@pytest.fixture(scope="module")
def t():
    return Target("q1:20:t1:20:-21.3:0.000944:1:       UUU           U:  UUUAG   GUUCGUUUAGA :  AAAUU   CAAGCAAAUCU :UA     U             C")


class TestTarget:
    def test_mfe(self, t):
        assert t.mfe == -21.3

    def test_pvalue(self, t):
        assert t.pvalue == 0.000944

    def test_position(self, t):
        assert t.position == 1

    def test_pairing(self, t):
        assert "\n".join(t.pairing) == """\
       UUU           U
  UUUAG   GUUCGUUUAGA 
  AAAUU   CAAGCAAAUCU 
UA     U             C"""

    def test_pairing_rc(self, t):
        assert "\n".join(t.pairing_rc) == """\
C             U     AU
 UCUAAACGAAC   UUAAA  
 AGAUUUGCUUG   GAUUU  
U           UUU       """

    def test_target(self, t):
        assert t.target == (
            'UUUAGUUUGUUCGUUUAGAU',
            '*****...***********.')

    def test_query(self, t):
        assert t.query == (
            'CUCUAAACGAACUUUAAAAU',
            '.***********.*****..')


@pytest.fixture(scope="module")
def t2():
    return Target("4:20:4:20:-12.7:0.934316:1:              UUU     U:       AUUUUAG   GUUCG :       UAAAAUU   CAAGC :GUGUGUC       U        ")


class TestTarget2:

    def test_position(self, t2):
        assert t2.position == 1

    def test_target(self, t2):
        assert t2.target == (
            'AUUUUAGUUUGUUCGU',
            '*******...*****.')

    def test_query(self, t2):
        assert t2.query == (
            'CGAACUUUAAAAUCUGUGUG',
            '*****.*******.......')

    def test_target_full(self, t2):
        t_full = "AUUUUAGUUUGUUCGUUUAG"
        assert t2.target_full(t_full) == (
            'AUUUUAGUUUGUUCGUUUAG',
            '*******...*****.....')

    def test_as_cofold(self, t2):
        t_full = "AUUUUAGUUUGUUCGUUUAG"
        assert t2.as_cofold(t_full) == (
            'CGAACUUUAAAAUCUGUGUG&AUUUUAGUUUGUUCGUUUAG',
            '(((((.(((((((.......&)))))))...))))).....')


@pytest.fixture(scope="module")
def t3():
    return Target('rgt:20:lft:20:-8.6:1.000000:4:     A    AAA    C     :      AUCU   ACAA      :      UAGA   UGUU      :CUUGUC           CUCUAG')


class TestTarget3:

    def test_position(self, t3):
        assert t3.position == 4

    def test_target(self, t3):
        assert t3.target == (
            'AAUCUAAAACAAC',
            '.****...****.')

    def test_query(self, t3):
        assert t3.query == (
            'GAUCUCUUGUAGAUCUGUUC',
            '......********......')

    def test_target_full(self, t3):
        t_full = "TGAAATCTAAAACAACACGA"
        assert t3.target_full(t_full) == (
            'UGAAAUCUAAAACAACACGA',
            '....****...****.....')

    def test_as_cofold(self, t3):
        t_full = "TGAAATCTAAAACAACACGA"
        assert t3.as_cofold(t_full) == (
            'GAUCUCUUGUAGAUCUGUUC&UGAAAUCUAAAACAACACGA',
            '......((((((((......&....))))...)))).....')


@pytest.fixture(scope="module")
def t4():
    return Target('t:20:q:20:-11.1:0.999999:1:            GUG     U :      UAGG C   ACAAG  :      GUCU G   UGUUC  :UCUCUU    A A       UC')


class TestTarget4:

    def test_position(self, t4):
        assert t4.position == 1

    def test_as_cofold(self, t4):
        t_full = "UAGGCGUGACAAGUUUCAUU"
        assert t4.as_cofold(t_full) == (
            'CUCUUGUAGAUCUGUUCUCU&UAGGCGUGACAAGUUUCAUU',
            '..(((((.(.((((......&)))))...))))).......')
