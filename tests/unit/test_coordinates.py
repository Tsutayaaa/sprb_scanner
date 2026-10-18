import pytest

from sprb_scanner.core.coordinates import (
    from_one_based_closed,
    overlap_length,
    to_one_based_closed,
)


@pytest.mark.parametrize("a,b", [(1, 1), (101, 200), (1, 100000)])
def test_coordinate_round_trip(a, b):
    interval = from_one_based_closed(a, b)
    assert interval == (a - 1, b)
    assert to_one_based_closed(*interval) == (a, b)


@pytest.mark.parametrize("a,b", [(0, 1), (2, 1), (True, 2), (1.5, 2)])
def test_invalid_external(a, b):
    with pytest.raises(ValueError):
        from_one_based_closed(a, b)


def test_overlap():
    assert overlap_length((0, 10), (10, 20)) == 0
    assert overlap_length((0, 10), (9, 20)) == 1
    with pytest.raises(ValueError):
        to_one_based_closed(10, 10)
