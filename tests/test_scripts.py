import subprocess as sp
from pathlib import Path

import pytest
from pygenomeviz.scripts.mmseqs import MMseqs
from pygenomeviz.scripts.mummer import MUMmer
from pygenomeviz.scripts.pmauve import ProgressiveMauve
from pygenomeviz.utils import load_dataset


def test_download_dataset_cli(tmp_path: Path):
    """Test download dataset CLI (pgv-download-dataset)"""
    dataset_name = "escherichia_phage"
    cmd = f"pgv-download-dataset -n {dataset_name} -o {tmp_path}"
    result = sp.run(cmd, shell=True)

    assert result.returncode == 0
    assert len(list(tmp_path.glob("*.gbk"))) == 4


def test_simpleplot_cli(tmp_path: Path):
    """Test simpleplot CLI (pgv-simpleplot)"""
    gbk_files, _ = load_dataset("enterobacteria_phage")
    gbk_files = [str(f) for f in gbk_files]

    # Test PNG,HTML image output
    for format in ("png", "html"):
        outfile = tmp_path / f"test.{format}"
        cmd = f"pgv-simpleplot --gbk_resources {' '.join((gbk_files))} -o {outfile}"
        result = sp.run(cmd, shell=True)
        assert result.returncode == 0
        assert outfile.exists()


@pytest.mark.skipif(
    condition=not MUMmer.check_installation(exit_on_false=False),
    reason="MUMmer is not installed in this environment.",
)
def test_mummer_cli(tmp_path: Path):
    """Test MUMmer CLI (pgv-mummer)"""
    gbk_files, _ = load_dataset("escherichia_phage")
    gbk_files = [str(f) for f in gbk_files]
    cmd = f"pgv-mummer --gbk_resources {' '.join(gbk_files)} -o {tmp_path} --curve"
    result = sp.run(cmd, shell=True)

    assert result.returncode == 0

    result_fig_file = tmp_path / "result.png"
    assert result_fig_file.exists()


@pytest.mark.skipif(
    condition=not MMseqs.check_installation(exit_on_false=False),
    reason="MMseqs is not installed in this environment.",
)
def test_mmseqs_cli(tmp_path: Path):
    """Test MMseqs CLI (pgv-mmseqs)"""
    gbk_files, _ = load_dataset("escherichia_phage")
    gbk_files = [str(f) for f in gbk_files]
    cmd = f"pgv-mmseqs --gbk_resources {' '.join(gbk_files)} -o {tmp_path} --curve"
    result = sp.run(cmd, shell=True)

    assert result.returncode == 0

    result_fig_file = tmp_path / "result.png"
    assert result_fig_file.exists()


@pytest.mark.skipif(
    condition=not ProgressiveMauve.check_installation(exit_on_false=False),
    reason="progressiveMauve is not installed in this environment.",
)
def test_pmauve_cli(tmp_path: Path):
    """Test progressiveMauve CLI (pgv-pmauve)"""
    gbk_files, _ = load_dataset("escherichia_phage")
    gbk_files = [str(f) for f in gbk_files]
    cmd = f"pgv-pmauve --seq_files {' '.join(gbk_files)} -o {tmp_path} --curve"
    result = sp.run(cmd, shell=True)

    assert result.returncode == 0

    result_fig_file = tmp_path / "result.png"
    assert result_fig_file.exists()
