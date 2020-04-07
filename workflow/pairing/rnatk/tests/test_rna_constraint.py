import pytest


from rnatk.rna.constraint import add_constraint, add_interconst, BPs


def test_bps():
    assert len(BPs) == 10


def test_add_constraint5():
    assert add_constraint('ACGGG&GGGGT', prime=5) == '((...&...))'
    assert add_constraint('ACGGG&GGGGT', prime=3) == '.....&.....'


def test_add_constraint3():
    assert add_constraint('GGGAC&GTGGG', prime=5) == '.....&.....'
    assert add_constraint('GGGAC&GTGGG', prime=3) == '...((&))...'


def test_add_constraint5ws():
    assert add_constraint('ACGGG&GGGGT', prime=5, struct='.(...&...).') == '((...&...))'
    assert add_constraint('ACGGG&GGGCG', prime=5, struct='.((..&...))') == '.((..&...))'


def test_add_constraint3ws():
    assert add_constraint('GGGAC&GTGGG', prime=3, struct='...(.&.)...') == '...((&))...'
    assert add_constraint('GGGAA&GTTGG', prime=3, struct='...((&.))..') == '...((&.))..'


def test_add_interconst():
    assert add_interconst("CCCCG&CGGGGT") == "<<<<<&>>>>>>"


def test_add_constraint3unbalance():
    assert add_constraint('CUCUUGUAGAUCUGUUCUCU&UAGGCGUGACAAGUUUCAUU',
                          prime=3,
                          struct='..(((((.(.((((......&)))))...))))).......') == \
           '..(((((.(.((((......&)))))...))))).......'


