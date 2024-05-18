from pathlib import Path

import pytest

from pygenomeviz.parser import Genbank


def test_genbank_property(gbk_file: Path):
    """Test Genbank instance properties"""
    gbk = Genbank(gbk_file)
    assert gbk.name == "test"
    assert len(gbk.records) == 1
    assert gbk.genome_length == 168903
    assert gbk.genome_length == gbk.full_genome_length
    assert gbk.genome_seq == gbk.full_genome_seq
    assert gbk.get_seqid2size() == {"NC_000866.4": 168903}
    assert gbk.extract_features() == list(gbk.get_seqid2features().values())[0]


def test_genbank_property_multi_record(multi_record_gbk_file: Path):
    """Test Genbank instance(multi-record) properties"""
    gbk = Genbank(multi_record_gbk_file)
    assert gbk.name == "multi_record"
    assert len(gbk.records) == 8
    assert gbk.genome_length == 136489
    assert gbk.full_genome_length == 1190241
    assert gbk.get_seqid2size() == {
        "NZ_LAEX01000001.1": 136489,
        "NZ_LAEX01000002.1": 155051,
        "NZ_LAEX01000003.1": 113854,
        "NZ_LAEX01000004.1": 41567,
        "NZ_LAEX01000005.1": 34901,
        "NZ_LAEX01000006.1": 132729,
        "NZ_LAEX01000007.1": 504069,
        "NZ_LAEX01000008.1": 71581,
    }
    assert gbk.extract_features() == list(gbk.get_seqid2features().values())[0]


def test_calc_genome_gc_content(multi_record_gbk_file: Path):
    """Test `calc_genome_gc_content()`"""
    gbk = Genbank(multi_record_gbk_file)
    default_gc = gbk.calc_genome_gc_content()
    assert default_gc == gbk.calc_genome_gc_content(seq=gbk.genome_seq)
    assert default_gc != gbk.calc_genome_gc_content(seq=gbk.full_genome_seq)

    seq = "ATGC" * 100
    assert gbk.calc_genome_gc_content(seq) == 50


def test_calc_gc_content(multi_record_gbk_file: Path):
    """Test `calc_gc_content()`"""
    gbk = Genbank(multi_record_gbk_file)
    default_gc_list = list(gbk.calc_gc_content()[1])
    assert default_gc_list == list(gbk.calc_gc_content(seq=gbk.genome_seq)[1])
    assert default_gc_list != list(gbk.calc_gc_content(seq=gbk.full_genome_seq)[1])

    assert max(gbk.calc_gc_content(seq="AT" * 10000)[1]) == 0
    assert max(gbk.calc_gc_content(seq="GC" * 10000)[1]) == 100


def test_calc_gc_skew(multi_record_gbk_file: Path):
    """Test `calc_gc_skew()`"""
    gbk = Genbank(multi_record_gbk_file)
    default_gc_skew_list = list(gbk.calc_gc_skew()[1])
    assert default_gc_skew_list == list(gbk.calc_gc_skew(seq=gbk.genome_seq)[1])
    assert default_gc_skew_list != list(gbk.calc_gc_skew(seq=gbk.full_genome_seq)[1])


def test_genbank_write_cds_fasta(multi_record_gbk_file: Path, tmp_path: Path):
    """Test write cds fasta result"""
    outfile = tmp_path / "cds.faa"
    Genbank(multi_record_gbk_file).write_cds_fasta(outfile)
    assert outfile.exists()


def test_genbank_write_genome_fasta(multi_record_gbk_file: Path, tmp_path: Path):
    """Test write genome fasta result"""
    # Check output file exists
    outfile = tmp_path / "genome.fna"
    Genbank(multi_record_gbk_file).write_genome_fasta(outfile)
    assert outfile.exists()

    # Check number of genome fasta
    with open(outfile) as f:
        lines = f.read().splitlines()
    fasta_count = len([line for line in lines if line.startswith(">")])
    assert fasta_count == 8


def test_parse_gbk_gz_file(gbk_gz_file: Path):
    """Parse GZ compressed genbank file"""
    gbk = Genbank(gbk_gz_file)
    assert gbk.name == "test"
    assert len(gbk.records) == 1


def test_parse_gbk_bz2_file(gbk_bz2_file: Path):
    """Parse BZ2 compressed genbank file"""
    gbk = Genbank(gbk_bz2_file)
    assert gbk.name == "test"
    assert len(gbk.records) == 1


def test_parse_zip_gbk_file(gbk_zip_file: Path):
    """Parse ZIP compressed genbank file"""
    gbk = Genbank(gbk_zip_file)
    assert gbk.name == "test"
    assert len(gbk.records) == 1


def test_parse_invalid_file_failed(gff_file: Path):
    """Test parse invalid(gff) file failed"""
    with pytest.raises(ValueError):
        Genbank(gff_file)
