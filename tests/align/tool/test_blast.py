from __future__ import annotations

from pathlib import Path

import pytest

from pygenomeviz.align import Blast
from tests.marker import skipif_blast_not_installed


def test_blast_get_tool_name():
    """Test `get_tool_name()`"""
    assert Blast.get_tool_name() == "BLAST"


def test_blast_get_binary_names():
    """Test `get_binary_names()`"""
    assert Blast.get_binary_names() == ["makeblastdb", "blastn", "tblastx"]


@skipif_blast_not_installed
@pytest.mark.parametrize(
    ("seqtype"),
    (
        pytest.param("nucleotide", id="blastn api (genbank dataset)"),
        pytest.param("protein", id="tblastx api (genbank dataset)"),
    ),
)
def test_blast_api_genbank(
    gbk_dataset_files: list[Path],
    tmp_path: Path,
    seqtype,
):
    """Run blast with genbank files"""
    align_coords = Blast(
        gbk_dataset_files,
        outdir=tmp_path,
        seqtype=seqtype,
    ).run()
    assert len(align_coords) > 0


@skipif_blast_not_installed
@pytest.mark.parametrize(
    ("seqtype"),
    (
        pytest.param("nucleotide", id="blastn api (fasta dataset)"),
        pytest.param("protein", id="tblastx api (fasta dataset)"),
    ),
)
def test_blast_api_fasta(
    fasta_dataset_files: list[Path],
    tmp_path: Path,
    seqtype,
):
    """Run blast with fasta files"""
    align_coords = Blast(
        fasta_dataset_files,
        outdir=tmp_path,
        seqtype=seqtype,
    ).run()
    assert len(align_coords) > 0
