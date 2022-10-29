from pathlib import Path

import pytest
from pygenomeviz import Gff
from pygenomeviz.parser.gff import GffRecord

###########################################################
# Test Gff
###########################################################


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


def _write_test_gff_file(
    test_gff_file: Path, delimiter: str = " ", annotation: bool = True
):
    gff_lines = (
        delimiter.join(("##sequence-region", "seqid1", "101", "1000")),
        "seqid1\tsource\tCDS\t501\t800\t.\t+\t.\tID=test1",
        delimiter.join(("##sequence-region", "seqid2", "201", "2000")),
        "seqid2\tsource\tCDS\t1001\t1300\t.\t+\t.\tID=test2",
        "seqid2\tsource\tCDS\t1401\t1820\t.\t+\t.\tID=test2",
    )
    if not annotation:
        gff_lines = [ln for ln in gff_lines if not ln.startswith("#")]
    with open(test_gff_file, "w") as f:
        f.write("\n".join(gff_lines))


def test_sequence_region_extract(tmp_path: Path):
    """Test sequence region extraction"""
    # Case1. Available sequence-region annotation line (delimiter: ' ' or '\t')
    gff_tmp_file = tmp_path / "test.gff"
    for delimiter in (" ", "\t"):
        _write_test_gff_file(gff_tmp_file, delimiter)
        gff = Gff(gff_tmp_file)
        assert gff.min_range == 100 and gff.max_range == 1000
        gff = Gff(gff_tmp_file, target_seqid="seqid2")
        assert gff.min_range == 200 and gff.max_range == 2000

    # Case2. Not available sequence-region annotation line
    _write_test_gff_file(gff_tmp_file, annotation=False)
    gff = Gff(gff_tmp_file)
    assert gff.min_range == 0 and gff.max_range == 800
    gff = Gff(gff_tmp_file, target_seqid="seqid2")
    assert gff.min_range == 0 and gff.max_range == 1820


def test_seqid_extract(tmp_path: Path):
    """Test seqid extraction"""
    gff_tmp_file = tmp_path / "test.gff"
    _write_test_gff_file(gff_tmp_file)

    # Case1. No specify target_seqid
    gff = Gff(gff_tmp_file)
    assert gff.target_seqid == "seqid1"
    assert gff.seqid_list == ["seqid1", "seqid2"]

    # Case2. Specify target_seqid
    gff = Gff(gff_tmp_file, target_seqid="seqid2")
    assert gff.target_seqid == "seqid2"
    assert gff.seqid_list == ["seqid1", "seqid2"]


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


###########################################################
# Test GffRecord
###########################################################


def test_to_gff_line():
    """Test GffRecord to_gff_line"""
    # Case1. score=None, strand=1, phase=None, attrs='ID=test'
    gff_record = GffRecord("1", "ncbi", "CDS", 10, 20, None, 1, None, {"ID": ["test"]})
    expected_gff_line = "1\tncbi\tCDS\t11\t20\t.\t+\t.\tID=test"
    assert gff_record.to_gff_line() == expected_gff_line

    # Case2. score=1e-10, strand=-1, phase=10, attrs='ID=test;Dbxref=ref1,ref2'
    attrs = {"ID": ["test"], "Dbxref": ["ref1", "ref2"]}
    gff_record = GffRecord("2", "ncbi", "mRNA", 20, 30, 1e-10, -1, 10, attrs)
    expected_gff_line = "2\tncbi\tmRNA\t21\t30\t1e-10\t-\t10\tID=test;Dbxref=ref1,ref2"
    assert gff_record.to_gff_line() == expected_gff_line


def test_parse_gff_line():
    """Test GffRecord parse_gff_line"""
    # Case1. score='.', strand='+', phase='.', attrs='ID=test'
    gff_line = "1\tncbi\tCDS\t11\t20\t.\t+\t.\tID=test"
    gff_record = GffRecord.parse_gff_line(gff_line)
    expected_gff_record = GffRecord(
        "1", "ncbi", "CDS", 10, 20, None, 1, None, {"ID": ["test"]}
    )
    assert gff_record == expected_gff_record

    # Case2. score=1e-10, strand='-', phase=10, attrs='ID=test;Dbxref=ref1,ref2'
    gff_line = "2\tncbi\tmRNA\t21\t30\t1e-10\t-\t10\tID=test;Dbxref=ref1,ref2"
    gff_record = GffRecord.parse_gff_line(gff_line)
    attrs = {"ID": ["test"], "Dbxref": ["ref1", "ref2"]}
    expected_gff_record = GffRecord("2", "ncbi", "mRNA", 20, 30, 1e-10, -1, 10, attrs)
    assert gff_record == expected_gff_record


def test_parse_and_to_gff_line():
    """Test `parse_gff_line` and `to_gff_line` results are compatible"""
    gff_line = "1\tncbi\tCDS\t11\t20\t.\t+\t.\tID=test"
    gff_record = GffRecord.parse_gff_line(gff_line)
    assert gff_line == gff_record.to_gff_line()

    gff_line = "2\tncbi\tmRNA\t21\t30\t1e-10\t-\t10\tID=test;Dbxref=ref1,ref2"
    gff_record = GffRecord.parse_gff_line(gff_line)
    assert gff_line == gff_record.to_gff_line()
