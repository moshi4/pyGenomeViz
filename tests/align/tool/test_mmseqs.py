from __future__ import annotations

from pathlib import Path

from pygenomeviz.align import MMseqs
from tests.marker import skipif_mmseqs_not_installed


def test_mmseqs_get_tool_name():
    """Test `get_tool_name()`"""
    assert MMseqs.get_tool_name() == "MMseqs"


def test_mmseqs_get_binary_names():
    """Test `get_binary_names()`"""
    assert MMseqs.get_binary_names() == ["mmseqs"]


@skipif_mmseqs_not_installed
def test_mmseqs_api(gbk_dataset_files: list[Path], tmp_path: Path):
    """Run mmseqs easy-rbh"""
    align_coords = MMseqs(gbk_dataset_files, outdir=tmp_path).run()
    assert len(align_coords) > 0
