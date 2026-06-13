from __future__ import annotations

import subprocess as sp
from typing import TYPE_CHECKING

from tests.marker import skipif_network_connection_failed

if TYPE_CHECKING:
    from pathlib import Path

CLI_NAME = "pgv-download"


@skipif_network_connection_failed
def test_download_cli(tmp_path: Path) -> None:
    """Run `pgv-download` cli"""
    cmd = f"{CLI_NAME} yersinia_phage -o {tmp_path}"
    cmd_res = sp.run(cmd, check=False, shell=True)
    assert cmd_res.returncode == 0
    assert len(list(tmp_path.glob("*.gbk"))) == 4


def test_download_cli_invalid_dataset_name_failed(tmp_path: Path) -> None:
    """Run `pgv-download` cli failed with invalid dataset name"""
    cmd = f"{CLI_NAME} invalid_dataset -o {tmp_path}"
    cmd_res = sp.run(cmd, check=False, shell=True)
    assert cmd_res.returncode == 2
