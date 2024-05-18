from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture
def testdata_dir() -> Path:
    """Testdata directory path"""
    return Path(__file__).parent / "testdata"


###########################################################
# Genbank test data
###########################################################


@pytest.fixture
def genbank_testdata_dir(testdata_dir: Path) -> Path:
    """Genbank testdata directory path"""
    return testdata_dir / "genbank"


@pytest.fixture
def gbk_file(genbank_testdata_dir: Path) -> Path:
    """Genbank file

    Escherichia phage T4: 1 record

    1. NC_000866.4 (168,903 bp)
    """
    return genbank_testdata_dir / "test.gbff"


@pytest.fixture
def gbk_gz_file(genbank_testdata_dir: Path) -> Path:
    """Genbank file (GZ compressed)"""
    return genbank_testdata_dir / "test.gbff.gz"


@pytest.fixture
def gbk_bz2_file(genbank_testdata_dir: Path) -> Path:
    """Genbank file (BZ2 compressed)"""
    return genbank_testdata_dir / "test.gbff.bz2"


@pytest.fixture
def gbk_zip_file(genbank_testdata_dir: Path) -> Path:
    """Genbank file (ZIP compressed)"""
    return genbank_testdata_dir / "test.zip"


@pytest.fixture
def multi_record_gbk_file(genbank_testdata_dir: Path) -> Path:
    """Multi-record genbank file

    Mycoplasma mycoides: 8 records

    1. NZ_LAEX01000001.1 (136,489 bp)
    2. NZ_LAEX01000002.1 (155,051 bp)
    2. NZ_LAEX01000003.1 (113,854 bp)
    3. NZ_LAEX01000004.1 ( 41,567 bp)
    4. NZ_LAEX01000005.1 ( 34,901 bp)
    5. NZ_LAEX01000006.1 (132,729 bp)
    6. NZ_LAEX01000007.1 (504,069 bp)
    7. NZ_LAEX01000008.1 ( 71,581 bp)
    """
    return genbank_testdata_dir / "multi_record.gbff"


@pytest.fixture
def gbk_dataset_files(genbank_testdata_dir: Path) -> list[Path]:
    """Genbank dataset files for genome comparison test

    influA viral genome dataset (Minimum dataset size)

    1. influA_California.gbff
    2. influA_Korea.gbff
    3. influA_NewYork.gbff
    4. influA_Shanghai.gbff
    """
    gbk_dataset_dir = genbank_testdata_dir / "influA"
    return list(gbk_dataset_dir.glob("*.gbff"))


###########################################################
# GFF test data
###########################################################


@pytest.fixture
def gff_testdata_dir(testdata_dir: Path) -> Path:
    """GFF testdata directory path"""
    return testdata_dir / "gff"


@pytest.fixture
def gff_file(gff_testdata_dir: Path) -> Path:
    """GFF file

    Escherichia phage T4: 1 record

    1. NC_000866.4 (168,903 bp)
    """
    return gff_testdata_dir / "test.gff"


@pytest.fixture
def gff_gz_file(gff_testdata_dir: Path) -> Path:
    """GFF file (GZ compressed)"""
    return gff_testdata_dir / "test.gff.gz"


@pytest.fixture
def gff_bz2_file(gff_testdata_dir: Path) -> Path:
    """GFF file (BZ2 compressed)"""
    return gff_testdata_dir / "test.gff.bz2"


@pytest.fixture
def gff_zip_file(gff_testdata_dir: Path) -> Path:
    """GFF file (ZIP compressed)"""
    return gff_testdata_dir / "test.zip"


@pytest.fixture
def multi_record_gff_file(gff_testdata_dir: Path) -> Path:
    """Multi-record GFF file

    Mycoplasma mycoides: 8 records

    1. NZ_LAEX01000001.1 (136,489 bp)
    2. NZ_LAEX01000002.1 (155,051 bp)
    2. NZ_LAEX01000003.1 (113,854 bp)
    3. NZ_LAEX01000004.1 ( 41,567 bp)
    4. NZ_LAEX01000005.1 ( 34,901 bp)
    5. NZ_LAEX01000006.1 (132,729 bp)
    6. NZ_LAEX01000007.1 (504,069 bp)
    7. NZ_LAEX01000008.1 ( 71,581 bp)
    """
    return gff_testdata_dir / "multi_record.gff"


###########################################################
# FASTA test data
###########################################################


@pytest.fixture
def fasta_testdata_dir(testdata_dir: Path) -> Path:
    """FASTA testdata directory path"""
    return testdata_dir / "fasta"


@pytest.fixture
def fasta_dataset_files(fasta_testdata_dir: Path) -> list[Path]:
    """FASTA dataset files for genome comparison test

    influA viral genome dataset (Minimum dataset size)

    1. influA_California.fna
    2. influA_Korea.fna
    3. influA_NewYork.fna
    4. influA_Shanghai.fna
    """
    fasta_dataset_dir = fasta_testdata_dir / "influA"
    return list(fasta_dataset_dir.glob("*.fna"))
