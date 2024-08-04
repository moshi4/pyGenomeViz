from __future__ import annotations

import logging
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Sequence

from pygenomeviz.align import AlignCoord
from pygenomeviz.align.tool import AlignToolBase
from pygenomeviz.parser import Fasta, Genbank


class ProgressiveMauve(AlignToolBase):
    """progressiveMauve Alignment Class"""

    def __init__(
        self,
        seqs: Sequence[str | Path | Fasta | Genbank],
        *,
        outdir: str | Path | None = None,
        refid: int = 0,
        cmd_opts: str | None = None,
        logger: logging.Logger | None = None,
        quiet: bool = True,
    ):
        """
        Parameters
        ----------
        seqs : Sequence[str | Path | Fasta | Genbank]
            List of fasta or genbank
            (file suffix must be `.fa`, `.fna`, `.fasta`, `.gb`, `.gbk`, `.gbff`)
        outdir : str | Path | None, optional
            Temporary result directory. If None, tmp directory is auto created.
        refid : int, optional
            Reference genome index
        cmd_opts : str | None, optional
            `progressiveMauve` additional command options
        logger : logging.Logger | None, optional
            Logger object. If None, logger instance newly created.
        quiet : bool, optional
            If True, don't display log message
        """
        super().__init__(logger, quiet)

        self._seqs = self._parse_input_gbk_and_fasta_seqs(seqs)
        self._outdir = None if outdir is None else Path(outdir)
        self._refid = refid
        self._cmd_opts = cmd_opts

    @classmethod
    def get_tool_name(cls) -> str:
        """Tool name"""
        return "progressiveMauve"

    @classmethod
    def get_binary_names(cls) -> list[str]:
        """Binary names"""
        return ["progressiveMauve"]

    @property
    def name2seqlen(self) -> dict[str, int]:
        """Name & sequence length dict"""
        name2seqlen = {}
        for seq in self._seqs:
            name2seqlen[seq.name] = sum(list(seq.get_seqid2size().values()))
        return name2seqlen

    def run(self) -> list[AlignCoord]:
        """Run genome alignment"""
        with TemporaryDirectory() as tmpdir:
            outdir = self._outdir if self._outdir else tmpdir
            outdir = Path(outdir)
            os.makedirs(outdir, exist_ok=True)
            genome_files: list[Path] = self._write_genome_files(outdir)

            # Run progressiveMauve
            xmfa_file = outdir / "pmauve.xmfa"
            bbone_file = outdir / "pmauve_bbone.tsv"
            self._logger.info(f"{'='*10} Start progressiveMauve Alignment {'='*10}")
            cmd = f"progressiveMauve --output={xmfa_file} --backbone-output={bbone_file} {' '.join(map(str, genome_files))}"  # noqa: E501
            if self._cmd_opts:
                cmd = f"{cmd} {self._cmd_opts}"
            self.run_cmd(cmd, self._logger)
            self._logger.info(f"{'='*10} Finish progressiveMauve Alignment {'='*10}")

            names = [file.stem for file in genome_files]
            return AlignCoord.parse_pmauve_file(bbone_file, names, self._refid)

    def _write_genome_files(self, outdir: str | Path) -> list[Path]:
        """Write genome fasta files to output directory

        Parameters
        ----------
        outdir : str | Path
            Target output directory

        Returns
        -------
        genome_files : list[Path]
            Genome fasta files
        """
        genome_files: list[Path] = []
        for seq in self._seqs:
            genome_file = Path(outdir) / f"{seq.name}.fna"
            cls_name = seq.__class__.__name__
            log_msg = f"Convert {cls_name} object to genome fasta file '{genome_file}'"
            self._logger.info(log_msg)
            seq.write_genome_fasta(genome_file)
            genome_files.append(genome_file)
        return genome_files
