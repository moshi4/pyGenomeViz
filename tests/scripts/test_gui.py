from __future__ import annotations

import subprocess as sp

CLI_NAME = "pgv-gui"


def test_gui_cli_help():
    """Run GUI cli with help option"""
    cmd = f"{CLI_NAME} --help"
    cmd_res = sp.run(cmd, shell=True)
    assert cmd_res.returncode == 0
