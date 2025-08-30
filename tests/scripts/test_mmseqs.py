from __future__ import annotations

import subprocess as sp
from typing import TYPE_CHECKING

from tests.marker import skipif_mmseqs_not_installed

if TYPE_CHECKING:
    from pathlib import Path

CLI_NAME = "pgv-mmseqs"


@skipif_mmseqs_not_installed
def test_mmseqs_cli(gbk_dataset_files: list[Path], tmp_path: Path) -> None:
    """Run CLI"""
    seqs = " ".join([str(file) for file in gbk_dataset_files])
    cmd = f"{CLI_NAME} {seqs} -o {tmp_path} --formats png"
    cmd_res = sp.run(cmd, check=False, shell=True)

    assert cmd_res.returncode == 0

    png_file = tmp_path / "result.png"
    assert png_file.exists()
