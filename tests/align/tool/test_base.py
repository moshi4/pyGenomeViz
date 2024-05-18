from __future__ import annotations

from pathlib import Path

import pytest

from pygenomeviz.align import Blast


def test_aligner_tool_initialize_failed(
    gbk_dataset_files: list[Path],
    monkeypatch,
):
    """Test aligner tool initialize failed when binary not found"""
    monkeypatch.setattr(Blast, "get_binary_names", lambda: ["pseudo_cmd"])
    with pytest.raises(RuntimeError):
        Blast(gbk_dataset_files)


def test_aligner_tool_check_installation_failed(monkeypatch):
    """Test aligner tool check installation failed when binary not found"""
    monkeypatch.setattr(Blast, "get_binary_names", lambda: ["pseudo_cmd"])
    with pytest.raises(RuntimeError):
        Blast.check_installation()
    assert Blast.check_installation(exit_on_false=False) is False
