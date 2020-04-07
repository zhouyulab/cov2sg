import pytest
from rnatk.util.rc import reverse_complement, reverse


def test_reverse_complement_1():
    assert reverse_complement('ACCTNAGGT') == "ACCTNAGGT"


def test_reverse_1():
    assert reverse("AAAACCCC") == "CCCCAAAA"

