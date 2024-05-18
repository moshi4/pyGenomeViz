from __future__ import annotations

from pathlib import Path

from pygenomeviz.align import ProgressiveMauve
from tests.marker import skipif_pmauve_not_installed


def test_pmauve_get_tool_name():
    """Test `get_tool_name()`"""
    assert ProgressiveMauve.get_tool_name() == "progressiveMauve"


def test_pmauve_get_binary_names():
    """Test `get_binary_names()`"""
    assert ProgressiveMauve.get_binary_names() == ["progressiveMauve"]


@skipif_pmauve_not_installed
def test_pmauve_api_genbank(gbk_dataset_files: list[Path], tmp_path: Path):
    """Run progressiveMauve with genbank files"""
    align_coords = ProgressiveMauve(gbk_dataset_files, outdir=tmp_path).run()
    assert len(align_coords) > 0


@skipif_pmauve_not_installed
def test_pmauve_api_fasta(fasta_dataset_files: list[Path], tmp_path: Path):
    """Run progressiveMauve with fasta files"""
    align_coords = ProgressiveMauve(fasta_dataset_files, outdir=tmp_path).run()
    assert len(align_coords) > 0
