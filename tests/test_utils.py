from pathlib import Path

import pytest
from pygenomeviz.utils import ColorCycler, load_dataset, load_example_gff


def test_load_dataset(tmp_path: Path):
    """Test load_dataset"""
    gbk_files, links = load_dataset("escherichia_phage", tmp_path)
    for gbk_file in gbk_files:
        assert gbk_file.exists()
    assert links != []


def test_load_dataset_invalid_name_error():
    """Test load_dataset no dataset error"""
    with pytest.raises(ValueError) as e:
        invalid_dataset_name = "nothing"
        load_dataset(invalid_dataset_name)
    assert "dataset not found" in str(e.value)


def test_load_example_gff(tmp_path: Path):
    """Test load_example_gff"""
    gff_file = load_example_gff("enterobacteria_phage.gff", tmp_path)
    assert gff_file.exists()


def test_load_example_gff_invalid_name_error():
    """Test load_example_gff no filename error"""
    with pytest.raises(ValueError) as e:
        invalid_filename = "nothing"
        load_example_gff(invalid_filename)
    assert "filename=" in str(e.value)


def test_color_cycler():
    """Test color cycler"""
    # Check get color list length
    assert len(ColorCycler.get_color_list()) == 10
    assert len(ColorCycler.get_color_list(5)) == 5
    assert len(ColorCycler.get_color_list(20)) == 20

    # Check cycle index, color
    assert ColorCycler(0) != ColorCycler(1)
    assert ColorCycler(0) == ColorCycler(10)
    assert ColorCycler(15) == ColorCycler(25)

    # Check cycle counter
    assert ColorCycler() != ColorCycler()
    assert ColorCycler.counter == 2

    # Check reset cycle
    ColorCycler.reset_cycle()
    assert ColorCycler.counter == 0

    # Check cmap change
    ColorCycler.set_cmap("tab20")
    with pytest.raises(ValueError):
        ColorCycler.set_cmap("invalid name")
    assert len(ColorCycler.get_color_list()) == 20
