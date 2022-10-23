import math
from pathlib import Path

import pytest
from Bio.Seq import reverse_complement
from pygenomeviz import Genbank


def test_default_param(gbk_file: Path):
    """Test default parameter result"""
    gbk = Genbank(gbk_file)
    assert gbk.name == "test"
    assert gbk.full_genome_length == gbk.genome_length == 43741
    assert gbk.full_genome_seq == gbk.genome_seq


def test_reverse_param(gbk_file: Path):
    """Test reverse parameter result"""
    normal_gbk = Genbank(gbk_file, reverse=False)
    reverse_gbk = Genbank(gbk_file, reverse=True)
    assert reverse_complement(normal_gbk.genome_seq) == reverse_gbk.genome_seq


def test_range_param(gbk_file: Path):
    """Test range parameter result"""
    min_range, max_range = 10000, 30000
    gbk = Genbank(gbk_file, min_range=min_range, max_range=max_range)
    expected_length = max_range - min_range

    assert gbk.full_genome_length != gbk.genome_length == expected_length
    assert gbk.full_genome_seq[min_range:max_range] == gbk.genome_seq


def test_range_error_param(gbk_file: Path):
    """Test range error parameter"""
    min_out_range, max_out_range = -100, 1000000
    # Case1. min_range < 0
    with pytest.raises(ValueError):
        Genbank(gbk_file, min_range=min_out_range, max_range=10000)
    # Case2. max_range > genome_length
    with pytest.raises(ValueError):
        Genbank(gbk_file, min_range=100, max_range=max_out_range)


def test_calc_genome_gc_content(gbk_file: Path):
    """Test genome GC content calculation"""
    gbk = Genbank(gbk_file)
    assert type(gbk.calc_genome_gc_content()) == float


def test_calc_gc_skew(gbk_file: Path):
    """Test GC skew calculation"""
    gbk = Genbank(gbk_file)
    # Default parameter
    pos_list, gc_skew_list = gbk.calc_gc_skew(None, None)
    expected_count = math.ceil(gbk.genome_length / int(gbk.genome_length / 1000)) + 1
    assert len(pos_list) == len(gc_skew_list) == expected_count
    # User setting parameter
    window_size, step_size = 500, 250
    pos_list, gc_skew_list = gbk.calc_gc_skew(window_size, step_size)
    expected_count = math.ceil(gbk.genome_length / step_size) + 1
    assert len(pos_list) == len(gc_skew_list) == expected_count


def test_calc_gc_content(gbk_file: Path):
    """Test GC content calculation"""
    gbk = Genbank(gbk_file)
    # Default parameter
    pos_list, gc_content_list = gbk.calc_gc_content(None, None)
    expected_count = math.ceil(gbk.genome_length / int(gbk.genome_length / 1000)) + 1
    assert len(pos_list) == len(gc_content_list) == expected_count
    # User setting parameter
    window_size, step_size = 500, 250
    pos_list, gc_content_list = gbk.calc_gc_content(window_size, step_size)
    expected_count = math.ceil(gbk.genome_length / step_size) + 1
    assert len(pos_list) == len(gc_content_list) == expected_count


def test_extract_features(gbk_file: Path):
    """Test write cds fasta"""
    gbk = Genbank(gbk_file)
    assert len(gbk.extract_features("CDS", None, False, False)) == 60


def test_write_cds_fasta(gbk_file: Path, tmp_path: Path):
    """Test write cds fasta"""
    cds_fasta_file = tmp_path / "cds.faa"
    gbk = Genbank(gbk_file)
    gbk.write_cds_fasta(cds_fasta_file)
    assert cds_fasta_file.exists()


def test_write_genome_fasta(gbk_file: Path, tmp_path: Path):
    """Test write genome fasta"""
    genome_fasta_file = tmp_path / "genome.fna"
    gbk = Genbank(gbk_file)
    gbk.write_genome_fasta(genome_fasta_file)
    assert genome_fasta_file.exists()


def test_parse_bzfile(gbk_bzfile: Path):
    """Test parse genbank file (bz2 compressed)"""
    gbk = Genbank(gbk_bzfile)
    assert gbk.name == "test"


def test_parse_gzfile(gbk_gzfile: Path):
    """Test parse genbank file (gz compressed)"""
    gbk = Genbank(gbk_gzfile)
    assert gbk.name == "test"


def test_parse_zipfile(gbk_zipfile: Path):
    """Test parse genbank file (zip compressed)"""
    gbk = Genbank(gbk_zipfile)
    assert gbk.name == "test"
