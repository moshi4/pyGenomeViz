from pathlib import Path

import pytest

###########################################################
# Genbank test data
###########################################################


@pytest.fixture
def gbk_testdata_dir() -> Path:
    """Genbank testdata directory"""
    return Path(__file__).parent / "testdata" / "genbank"


@pytest.fixture
def gbk_file(gbk_testdata_dir: Path) -> Path:
    """Genbank file"""
    return gbk_testdata_dir / "test.gbff"


@pytest.fixture
def gbk_bzfile(gbk_testdata_dir: Path) -> Path:
    """Genbank file (bz2 compressed)"""
    return gbk_testdata_dir / "test.gbff.bz2"


@pytest.fixture
def gbk_gzfile(gbk_testdata_dir: Path) -> Path:
    """Genbank file (gz compressed)"""
    return gbk_testdata_dir / "test.gbff.gz"


@pytest.fixture
def gbk_zipfile(gbk_testdata_dir: Path) -> Path:
    """Genbank file (zip compressed)"""
    return gbk_testdata_dir / "test.zip"


###########################################################
# GFF test data
###########################################################


@pytest.fixture
def gff_testdata_dir() -> Path:
    """GFF testdata directory"""
    return Path(__file__).parent / "testdata" / "gff"


@pytest.fixture
def gff_file(gff_testdata_dir: Path) -> Path:
    """GFF file"""
    return gff_testdata_dir / "test.gff"


@pytest.fixture
def gff_bzfile(gff_testdata_dir: Path) -> Path:
    """GFF file (bz2 compressed)"""
    return gff_testdata_dir / "test.gff.bz2"


@pytest.fixture
def gff_gzfile(gff_testdata_dir: Path) -> Path:
    """GFF file (gz compressed)"""
    return gff_testdata_dir / "test.gff.gz"


@pytest.fixture
def gff_zipfile(gff_testdata_dir: Path) -> Path:
    """GFF file (zip compressed)"""
    return gff_testdata_dir / "test.zip"
