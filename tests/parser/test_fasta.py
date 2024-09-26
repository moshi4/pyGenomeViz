from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from pygenomeviz.parser import Fasta

fasta_txt = textwrap.dedent(
    """
    >species1
    A
    >species2_1
    TT
    >species2_2
    CCC
    >species3 annotation
    GGGG
    """
)[1:]


def test_fasta_property(tmp_path: Path):
    """Test Fasta instance properties"""
    fasta_file = tmp_path / "test.fa"
    with open(fasta_file, "w", encoding="utf-8") as f:
        f.write(fasta_txt)
    fasta = Fasta(fasta_file)

    assert fasta.name == "test"
    assert fasta.genome_seq == "A"
    assert fasta.genome_length == 1
    assert fasta.full_genome_seq == "ATTCCCGGGG"
    assert fasta.full_genome_length == 10

    assert fasta.get_seqid2seq() == dict(
        species1="A",
        species2_1="TT",
        species2_2="CCC",
        species3="GGGG",
    )
    assert fasta.get_seqid2size() == dict(
        species1=1,
        species2_1=2,
        species2_2=3,
        species3=4,
    )


def test_parse_fasta_gz_file(fasta_gz_file: Path):
    """Parse GZ compressed fasta file"""
    fasta = Fasta(fasta_gz_file)
    assert fasta.name == "test"
    assert len(fasta.records) == 8


def test_parse_fasta_bz2_file(fasta_bz2_file: Path):
    """Parse bz2 compressed fasta file"""
    fasta = Fasta(fasta_bz2_file)
    assert fasta.name == "test"
    assert len(fasta.records) == 8


def test_parse_fasta_zip_file(fasta_zip_file: Path):
    """Parse zip compressed fasta file"""
    fasta = Fasta(fasta_zip_file)
    assert fasta.name == "test"
    assert len(fasta.records) == 8


def test_parse_invalid_file_failed(gff_file: Path):
    """Test parse invalid file failed"""
    with pytest.raises(ValueError):
        Fasta(gff_file)
