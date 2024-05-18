from __future__ import annotations

import subprocess as sp
from pathlib import Path

from tests.marker import skipif_mmseqs_not_installed

CLI_NAME = "pgv-mmseqs"


@skipif_mmseqs_not_installed
def test_mmseqs_cli(gbk_dataset_files: list[Path], tmp_path: Path):
    """Run CLI"""
    seqs = " ".join([str(file) for file in gbk_dataset_files])
    cmd = f"{CLI_NAME} {seqs} -o {tmp_path} --formats png"
    cmd_res = sp.run(cmd, shell=True)

    assert cmd_res.returncode == 0

    png_file = tmp_path / "result.png"
    assert png_file.exists()
