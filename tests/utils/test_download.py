from pathlib import Path

import pytest

from pygenomeviz.utils import (
    fetch_genbank_by_accid,
    load_example_genbank_dataset,
    load_example_gff_file,
)
from tests.marker import skipif_network_connection_failed


@skipif_network_connection_failed
def test_load_example_genbank_dataset_success():
    """Test `load_example_genbank_dataset()`"""
    gbk_files = load_example_genbank_dataset("enterobacteria_phage")
    assert all([gbk_file.exists() for gbk_file in gbk_files])


def test_load_example_genbank_dataset_fail():
    """Test `load_example_genbank_dataset()` with invalid dataset name"""
    with pytest.raises(ValueError):
        load_example_genbank_dataset("invalid_name")  # type: ignore


@skipif_network_connection_failed
def test_load_example_gff_file_success():
    """Test `load_example_gff_file()`"""
    gff_file = load_example_gff_file("enterobacteria_phage.gff")
    assert gff_file.exists()


def test_load_example_gff_file_fail():
    """Test `load_example_gff_file()` with invalid name"""
    with pytest.raises(ValueError):
        load_example_gff_file("invalid_name")  # type: ignore


@pytest.mark.skip("Requires network connection and unstable test result")
def test_fetch_genbank_by_accid(tmp_path: Path):
    """Test `fetch_genbank_by_accid()`"""
    outfile = tmp_path / "download.gbk"
    fetch_genbank_by_accid("NC_013600", gbk_outfile=outfile)
    assert outfile.exists()
