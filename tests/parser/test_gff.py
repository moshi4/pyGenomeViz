from pathlib import Path

import pytest
from pygenomeviz import Gff


def test_gff_load_default(gff_file: Path):
    """Test gff load default"""
    gff = Gff(gff_file)
    assert gff.name == "test"
    assert gff.min_range == 0
    assert gff.max_range == 60942
    assert gff.range_size == 60942


def test_gff_load_user_params(gff_file: Path):
    """Test gff load user params"""
    name, seqid, min_range, max_range = "user", "NC_000902.1", 10000, 30000
    gff = Gff(
        gff_file,
        name=name,
        target_seqid=seqid,
        min_range=min_range,
        max_range=max_range,
    )
    assert gff.name == name
    assert gff.min_range == min_range
    assert gff.max_range == max_range
    assert gff.range_size == max_range - min_range
    assert seqid == gff.seqid_list[0]


def test_gff_load_seqid_exist_error(gff_file: Path):
    """Test gff load seqid exist error"""
    with pytest.raises(ValueError):
        Gff(gff_file, target_seqid="no_exists")


def test_gff_load_range_error(gff_file: Path):
    """Test gff load range error"""
    # Case1. min_range < 0
    with pytest.raises(ValueError):
        Gff(gff_file, min_range=-10000, max_range=40000)
    # Case2. min_range > max_range
    with pytest.raises(ValueError):
        Gff(gff_file, min_range=40000, max_range=10000)


def test_extract_features(gff_file: Path):
    """Test extract_features"""
    gff = Gff(gff_file)
    cds_features = gff.extract_features()
    assert len(cds_features) == 83

    no_features = gff.extract_features("no_data")
    assert len(no_features) == 0


def test_extract_features_restrict_range(gff_file: Path):
    """Test extract_features with restrict range"""
    gff = Gff(gff_file, min_range=10000, max_range=30000)
    features = gff.extract_features()
    assert len(features) == 28


def test_parse_bzfile(gff_bzfile: Path):
    """Test parse GFF file (bz2 compressed)"""
    gff = Gff(gff_bzfile)
    assert gff.name == "test"


def test_parse_gzfile(gff_gzfile: Path):
    """Test parse GFF file (gz compressed)"""
    gff = Gff(gff_gzfile)
    assert gff.name == "test"


def test_parse_zipfile(gff_zipfile: Path):
    """Test parse GFF file (zip compressed)"""
    gff = Gff(gff_zipfile)
    assert gff.name == "test"
