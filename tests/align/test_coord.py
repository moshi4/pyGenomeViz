from __future__ import annotations

from pathlib import Path

import pytest

from pygenomeviz.align import AlignCoord


@pytest.mark.parametrize(
    ("qstart", "qend", "rstart", "rend", "qstrand", "rstrand", "is_inverted"),
    (
        pytest.param(0, 100, 100, 200, 1, 1, False, id="forward-forward align coord"),
        pytest.param(0, 100, 200, 100, 1, -1, True, id="forward-reverse align coord"),
        pytest.param(100, 0, 200, 100, -1, -1, False, id="reverse-reverse align coord"),
    ),
)
def test_align_coord_property(
    qstart: int,
    qend: int,
    rstart: int,
    rend: int,
    qstrand: int,
    rstrand: int,
    is_inverted: bool,
):
    """Test `AlignCoord` property"""
    qid, qname, qlen = "qid", "qname", abs(qend - qstart)
    rid, rname, rlen = "rid", "rname", abs(rend - rstart)
    identity, evalue = 95.5, 1e-20

    ac = AlignCoord(
        qid,
        qname,
        qstart,
        qend,
        rid,
        rname,
        rstart,
        rend,
        identity,
        evalue,
    )

    assert ac.query_id == qid
    assert ac.query_name == qname
    assert ac.query_start == qstart
    assert ac.query_end == qend
    assert ac.query_length == qlen
    assert ac.ref_id == rid
    assert ac.ref_name == rname
    assert ac.ref_start == rstart
    assert ac.ref_end == rend
    assert ac.ref_length == rlen
    assert ac.identity == identity
    assert ac.evalue == evalue
    assert ac.query_strand == qstrand
    assert ac.ref_strand == rstrand
    assert ac.is_inverted == is_inverted

    tsv_format_line = f"{qid}\t{qname}\t{qstart}\t{qend}\t{qlen}\t{rid}\t{rname}\t{rstart}\t{rend}\t{rlen}\t{identity}\t{evalue}"  # noqa: E501
    assert ac.as_tsv_format == tsv_format_line


@pytest.fixture()
def align_coords() -> list[AlignCoord]:
    """List of AlignCoord"""
    align_coords = []
    ac1 = AlignCoord("q", "q1", 0, 10, "r", "r1", 5, 15, None, None)
    ac2 = AlignCoord("q", "q2", 0, 20, "r", "r2", 5, 25, 80.0, None)
    ac3 = AlignCoord("q", "q3", 0, 30, "r", "r3", 35, 5, None, 1e-20)
    ac4 = AlignCoord("q", "q4", 0, 40, "r", "r4", 45, 5, 50.0, 1e-3)
    align_coords.extend((ac1, ac2, ac3, ac4))
    return align_coords


def test_align_coord_write_read(align_coords: list[AlignCoord], tmp_path: Path):
    """Test `AlignCoord.write()` & `AlignCoord.read()` result"""
    # Test write
    outfile = tmp_path / "align_coords.tsv"
    AlignCoord.write(align_coords, outfile)
    assert outfile.exists()

    # Test read
    read_align_coords = AlignCoord.read(outfile)

    # Test write & read result is same
    for ac, read_ac in zip(align_coords, read_align_coords):
        assert ac == read_ac


@pytest.mark.parametrize(
    ("length_thr", "identity_thr", "evalue_thr", "expected_count"),
    (
        pytest.param(None, None, None, 4, id="No filter params"),
        pytest.param(20, None, None, 3, id="Length filter case1"),
        pytest.param(50, None, None, 0, id="Length filter case2"),
        pytest.param(None, 80, None, 3, id="Identity filter"),
        pytest.param(None, None, 1e-100, 2, id="E-value filter"),
        pytest.param(20, 90, 1e-10, 1, id="All filter params"),
    ),
)
def test_align_coord_filter(
    align_coords: list[AlignCoord],
    length_thr: int | None,
    identity_thr: float | None,
    evalue_thr: float | None,
    expected_count: int,
):
    """Test `AlignCoord.filter()` count"""
    filter_align_coords = AlignCoord.filter(
        align_coords,
        length_thr=length_thr,
        identity_thr=identity_thr,
        evalue_thr=evalue_thr,
    )
    assert len(filter_align_coords) == expected_count


def test_contains_special_method():
    """Test `__contains__` overlap check special method"""
    # Case1. Same object
    ac1 = AlignCoord("q1", "q1", 10, 20, "r1", "r1", 10, 20)
    assert (ac1 in ac1) is True
    # Case2. Both query, ref within
    ac2 = AlignCoord("q1", "q1", 0, 30, "r1", "r1", 0, 30)
    assert (ac1 in ac2) is True
    # Case3. Only ref within
    ac3 = AlignCoord("q1", "q1", 15, 30, "r1", "r1", 0, 30)
    assert (ac1 in ac3) is False
    # Case4. Only query within
    ac4 = AlignCoord("q1", "q1", 0, 30, "r1", "r1", 0, 15)
    assert (ac1 in ac4) is False
    # Case5. No overlap
    ac5 = AlignCoord("q1", "q1", 100, 50, "r1", "r1", 80, 30)
    assert (ac1 in ac5) is False
