from pathlib import Path

import pytest


@pytest.fixture
def testdata_dir() -> Path:
    """Testdata directory"""
    return Path(__file__).parent / "testdata"


@pytest.fixture
def gbk_file(testdata_dir: Path) -> Path:
    """Genbank file"""
    return testdata_dir / "test.gbff"


@pytest.fixture
def gbk_bzfile(testdata_dir: Path) -> Path:
    """Genbank file (bz2 compressed)"""
    return testdata_dir / "test.gbff.bz2"


@pytest.fixture
def gbk_gzfile(testdata_dir: Path) -> Path:
    """Genbank file (gz compressed)"""
    return testdata_dir / "test.gbff.gz"


@pytest.fixture
def gbk_zipfile(testdata_dir: Path) -> Path:
    """Genbank file (zip compressed)"""
    return testdata_dir / "test.zip"
