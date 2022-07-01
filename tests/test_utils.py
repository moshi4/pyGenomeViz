from pathlib import Path

import pytest
from pygenomeviz.utils import ColorCycler, load_dataset


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


def test_color_cycler():
    """Test color cycler"""
    assert ColorCycler(0) != ColorCycler(1)
    assert ColorCycler(0) == ColorCycler(10)
    assert ColorCycler(15) == ColorCycler(25)

    assert ColorCycler() != ColorCycler()
    assert ColorCycler.counter == 2

    ColorCycler.reset_cycle()
    assert ColorCycler.counter == 0

    ColorCycler.set_cmap("tab20")
    with pytest.raises(ValueError):
        ColorCycler.set_cmap("invalid name")
