"""Conversions at input/output boundaries; internal intervals are [start, end)."""


def validate_interval(start: int, end: int) -> None:
    if type(start) is not int or type(end) is not int:
        raise ValueError("Coordinates must be integers")
    if start < 0 or end <= start:
        raise ValueError("Expected a nonempty 0-based half-open interval")


def from_one_based_closed(start: int, end: int) -> tuple[int, int]:
    if type(start) is not int or type(end) is not int or start < 1 or end < start:
        raise ValueError("Expected a nonempty 1-based closed interval")
    return start - 1, end


def to_one_based_closed(start: int, end: int) -> tuple[int, int]:
    validate_interval(start, end)
    return start + 1, end


def overlap_length(a: tuple[int, int], b: tuple[int, int]) -> int:
    validate_interval(*a)
    validate_interval(*b)
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))
