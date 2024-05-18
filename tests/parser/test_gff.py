from pathlib import Path

import pytest

from pygenomeviz.parser import Gff


def test_gff_property(gff_file: Path):
    """Test Gff instance properties"""
    gff = Gff(gff_file)
    assert gff.name == "test"
    assert gff.seq_region == (0, 168903)
    assert gff.genome_length == 168903
    assert gff.genome_length == gff.full_genome_length
    assert gff.target_seqid == "NC_000866.4"
    assert gff.seqid_list == ["NC_000866.4"]
    assert len(gff.records) == 604
    assert len(gff.records) == len(gff.all_records)
    assert gff.get_seqid2size() == {"NC_000866.4": 168903}
    assert gff.extract_features() == list(gff.get_seqid2features().values())[0]


def test_gff_property_multi_record(multi_record_gff_file: Path):
    """Test Gff instance(multi-record) properties"""
    gff = Gff(multi_record_gff_file)
    assert gff.name == "multi_record"
    assert gff.seq_region == (0, 136489)
    assert gff.genome_length == 136489
    assert gff.full_genome_length == 1190241
    assert len(gff.records) == 246
    assert len(gff.all_records) == 2255

    expected_seqid2size = {
        "NZ_LAEX01000001.1": 136489,
        "NZ_LAEX01000002.1": 155051,
        "NZ_LAEX01000003.1": 113854,
        "NZ_LAEX01000004.1": 41567,
        "NZ_LAEX01000005.1": 34901,
        "NZ_LAEX01000006.1": 132729,
        "NZ_LAEX01000007.1": 504069,
        "NZ_LAEX01000008.1": 71581,
    }
    assert gff.target_seqid == list(expected_seqid2size.keys())[0]
    assert gff.seqid_list == list(expected_seqid2size.keys())
    assert gff.get_seqid2size() == expected_seqid2size
    assert gff.extract_features() == list(gff.get_seqid2features().values())[0]


def test_parse_gff_gz_file(gff_gz_file: Path):
    """Parse GZ compressed gff file"""
    gff = Gff(gff_gz_file)
    assert gff.name == "test"
    assert len(gff.records) == 604


def test_parse_gff_bz2_file(gff_bz2_file: Path):
    """Parse BZ2 compressed gff file"""
    gff = Gff(gff_bz2_file)
    assert gff.name == "test"
    assert len(gff.records) == 604


def test_parse_gff_zip_file(gff_zip_file: Path):
    """Parse ZIP compressed gff file"""
    gff = Gff(gff_zip_file)
    assert gff.name == "test"
    assert len(gff.records) == 604


def test_parse_invalid_file_failed(gbk_file: Path):
    """Test parse invalid(genbank) file failed"""
    with pytest.raises(ValueError):
        Gff(gbk_file)
