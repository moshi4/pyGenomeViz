from __future__ import annotations

import subprocess as sp
from pathlib import Path

import pytest

from tests.marker import skipif_mummer_not_installed

CLI_NAME = "pgv-mummer"


@skipif_mummer_not_installed
@pytest.mark.parametrize(
    ("seqtype"),
    (
        pytest.param("nucleotide", id="nucmer cli"),
        pytest.param("protein", id="promer cli"),
    ),
)
def test_mummer_cli(
    gbk_dataset_files: list[Path],
    tmp_path: Path,
    seqtype,
):
    """Run mummer cli"""
    seqs = " ".join([str(file) for file in gbk_dataset_files])
    cmd = f"{CLI_NAME} {seqs} -o {tmp_path} --formats png --seqtype {seqtype}"
    cmd_res = sp.run(cmd, shell=True)

    assert cmd_res.returncode == 0

    png_file = tmp_path / "result.png"
    assert png_file.exists()
