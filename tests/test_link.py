import pytest
from pygenomeviz.link import Link


def test_is_inverted():
    """Test inverted check"""
    # Case1. forward-forward
    ff_link = Link("1", 1, 100, "2", 101, 200)
    assert ff_link.is_inverted is False
    # Case2. forward-reverse
    fr_link = Link("1", 1, 100, "2", 200, 101)
    assert fr_link.is_inverted is True
    # Case3. reverse-reverse
    rr_link = Link("1", 100, 1, "2", 200, 101)
    assert rr_link.is_inverted is False


def test_add_offset():
    """Test add offset to link"""
    name1, start1, end1, offset1 = "link1", 1, 100, 200
    name2, start2, end2, offset2 = "link2", 101, 200, 250
    link = Link(name1, start1, end1, name2, start2, end2)
    offset_link = link.add_offset({name1: offset1, name2: offset2})

    # Check original instance is not changed
    assert link.track_start1 == start1
    assert link.track_end1 == end1
    assert link.track_start2 == start2
    assert link.track_end2 == end2
    # Check offset added new instance is changed
    assert offset_link.track_start1 == start1 + offset1
    assert offset_link.track_end1 == end1 + offset1
    assert offset_link.track_start2 == start2 + offset2
    assert offset_link.track_end2 == end2 + offset2


def test_color_string_error():
    """Test color string error"""
    link1, link2 = ("link1", 1, 100), ("link2", 101, 200)
    invalid_color = "nocolor"
    # Case1. normal color is not color like string
    with pytest.raises(ValueError):
        Link(*link1, *link2, normal_color=invalid_color)
    # Case2. inverted color is not color like string
    with pytest.raises(ValueError):
        Link(*link1, *link2, inverted_color=invalid_color)


def test_size_ratio_error():
    """Test size ratio error"""
    link1, link2 = ("link1", 1, 100), ("link2", 101, 200)
    # Case1. size_ratio < 0
    with pytest.raises(ValueError):
        Link(*link1, *link2, size_ratio=-1)
    # Case2. size_ratio > 1
    with pytest.raises(ValueError):
        Link(*link1, *link2, size_ratio=2)


def test_interpolation_value_error():
    """Test interpolation value error"""
    link1, link2 = ("link1", 1, 100), ("link2", 101, 200)
    # Case1. vmin < 0
    with pytest.raises(ValueError):
        Link(*link1, *link2, v=50, vmin=-100)
    # Case2. vmax > 100
    with pytest.raises(ValueError):
        Link(*link1, *link2, v=50, vmax=200)
    # Case3. v < vmin
    with pytest.raises(ValueError):
        Link(*link1, *link2, v=50, vmin=60)
    # Case4. v > vmax
    with pytest.raises(ValueError):
        Link(*link1, *link2, v=90, vmax=80)
