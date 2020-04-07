import pytest
from rnatk.ios.fastq import reader_fastq
from rnatk.util.rc import reverse_complement
from rnatk.op.pemerge import semiglobal, align_and_correct


@pytest.fixture()
def pereads():
    fq1 = "tests/data/slamtrim_R1.fq"
    fq2 = "tests/data/slamtrim_R2.fq"
    yield zip(reader_fastq(fq1), reader_fastq(fq2))


def test_semiglobal(pereads):
    pes = [(r1.get_seq(), reverse_complement(r2.get_seq())) for (r1, r2) in pereads]
    res = semiglobal(*pes[0])
    assert res[2] == 58 and res[4] == 75


def test_align_and_correct(pereads):
    pes = [(r1, r2) for (r1, r2) in pereads]
    s1, q1 = pes[0][0].get_seq(), pes[0][0].get_qual()
    s2, q2 = reverse_complement(pes[0][1].get_seq()), "".join(reversed(pes[0][1].get_qual()))
    s1c, q1c, s2c, q2c, match_ratio = align_and_correct(s1, q1, s2, q2)
    assert match_ratio > 0.89 and s1c[1] == s2[0] and q1c[1] == q2[0]
