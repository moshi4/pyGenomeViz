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
    link_name1, link_name2 = "link1", "link2"
    link = Link(link_name1, 1, 100, link_name2, 101, 200)
    offset_link = link.add_offset({link_name1: 200, link_name2: 250})

    # Check original instance is not changed
    assert link.track_start1 == 1 and link.track_end1 == 100
    assert link.track_start2 == 101 and link.track_end2 == 200
    # Check offset added new instance is changed
    assert offset_link.track_start1 == 201 and offset_link.track_end1 == 300
    assert offset_link.track_start2 == 351 and offset_link.track_end2 == 450
