from __future__ import annotations

import subprocess as sp
from pathlib import Path

import pytest

from tests.marker import skipif_blast_not_installed

CLI_NAME = "pgv-blast"


@skipif_blast_not_installed
@pytest.mark.parametrize(
    ("seqtype"),
    (
        pytest.param("nucleotide", id="blastn cli"),
        pytest.param("protein", id="tblastx cli"),
    ),
)
def test_blast_cli(
    gbk_dataset_files: list[Path],
    tmp_path: Path,
    seqtype,
):
    """Test blast cli"""
    seqs = " ".join([str(file) for file in gbk_dataset_files])
    cmd = f"{CLI_NAME} {seqs} -o {tmp_path} --formats png --seqtype {seqtype}"
    cmd_res = sp.run(cmd, shell=True)

    assert cmd_res.returncode == 0

    png_file = tmp_path / "result.png"
    assert png_file.exists()


@pytest.mark.parametrize(
    ["cmd"],
    [
        pytest.param(
            f"{CLI_NAME} seq1.gbk",
            id="Less than two files input (seqs positional args)",
        ),
        pytest.param(
            f"{CLI_NAME} seq1.gbk seq2.gbk --normal_link_color invalid_color",
            id="Invalid color name is set in --normal_link_color option",
        ),
        pytest.param(
            f"{CLI_NAME} seq1.gbk seq2.gbk --inverted_link_color invalid_color",
            id="Invalid color name is set in --inverted_link_color option",
        ),
        pytest.param(
            f"{CLI_NAME} seq1.gbk seq2.gbk --feature_type2color invalid",
            id="Invalid argument style of --feature_type2color option",
        ),
        pytest.param(
            f"{CLI_NAME} seq1.gbk seq2.gbk --feature_type2color CDS:orange tRNA:test",
            id="Invalid color name is set in --feature_type2color option",
        ),
    ],
)
def test_cli_parser_error_case(cmd: str, tmp_path: Path):
    """Test cli parser error cases"""
    cmd += f" -o {tmp_path}"
    cmd_res = sp.run(cmd, shell=True)
    assert cmd_res.returncode == 2
