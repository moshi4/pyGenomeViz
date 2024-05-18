from __future__ import annotations

from pathlib import Path

import pytest

from pygenomeviz.align import MUMmer
from tests.marker import skipif_mummer_not_installed


def test_mummer_get_tool_name():
    """Test `get_tool_name()`"""
    assert MUMmer.get_tool_name() == "MUMmer"


def test_mummer_get_binary_names():
    """Test `get_binary_names()`"""
    assert MUMmer.get_binary_names() == [
        "nucmer",
        "promer",
        "delta-filter",
        "show-coords",
    ]


@skipif_mummer_not_installed
@pytest.mark.parametrize(
    ("seqtype"),
    (
        pytest.param("nucleotide", id="nucmer api (genbank dataset)"),
        pytest.param("protein", id="promer api (genbank dataset)"),
    ),
)
def test_mummer_api_genbank(
    gbk_dataset_files: list[Path],
    tmp_path: Path,
    seqtype,
):
    """Run mummer with genbank files"""
    align_coords = MUMmer(
        gbk_dataset_files,
        outdir=tmp_path,
        seqtype=seqtype,
    ).run()
    assert len(align_coords) > 0


@skipif_mummer_not_installed
@pytest.mark.parametrize(
    ("seqtype"),
    (
        pytest.param("nucleotide", id="nucmer api (fasta dataset)"),
        pytest.param("protein", id="promer api (fasta dataset)"),
    ),
)
def test_mummer_api_fasta(
    fasta_dataset_files: list[Path],
    tmp_path: Path,
    seqtype,
):
    """Run mummer with fasta files"""
    align_coords = MUMmer(
        fasta_dataset_files,
        outdir=tmp_path,
        seqtype=seqtype,
    ).run()
    assert len(align_coords) > 0
