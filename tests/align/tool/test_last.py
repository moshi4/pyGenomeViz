from __future__ import annotations

from pathlib import Path

from pygenomeviz.align import Last
from pygenomeviz.const import UNKNOWN_VERSION
from tests.marker import skipif_last_not_installed


def test_last_get_tool_name():
    """Test `get_tool_name()`"""
    assert Last.get_tool_name() == "Last"


def test_last_get_binary_name():
    """Test `get_binary_names()`"""
    assert Last.get_binary_names() == ["lastdb", "lastal", "last-split", "maf-convert"]


@skipif_last_not_installed
def test_last_get_version():
    """Test `get_version()`"""
    assert Last.get_version() != UNKNOWN_VERSION


@skipif_last_not_installed
def test_last_api_genbank(
    gbk_dataset_files: list[Path],
    tmp_path: Path,
):
    """Run Last with genbank files"""
    align_coords = Last(
        gbk_dataset_files,
        outdir=tmp_path,
    ).run()
    assert len(align_coords) > 0


@skipif_last_not_installed
def test_last_api_fasta(
    fasta_dataset_files: list[Path],
    tmp_path: Path,
):
    """Run Last with fasta files"""
    align_coords = Last(
        fasta_dataset_files,
        outdir=tmp_path,
    ).run()
    assert len(align_coords) > 0
