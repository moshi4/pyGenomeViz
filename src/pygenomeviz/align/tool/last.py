from __future__ import annotations

import logging
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Sequence

from pygenomeviz.align import AlignCoord
from pygenomeviz.align.tool import AlignToolBase
from pygenomeviz.parser import Fasta, Genbank

logger = logging.getLogger(__name__)


class Last(AlignToolBase):
    """Last Alignment Class

    This class is experimental. API may change in the future release.
    """

    def __init__(
        self,
        seqs: Sequence[str | Path | Fasta | Genbank],
        *,
        outdir: str | Path | None = None,
        threads: int | None = None,
        cmd_opts: str | None = None,
    ):
        """
        Parameters
        ----------
        seqs : Sequence[str | Path | Fasta | Genbank]
            List of fasta or genbank
            (file suffix must be `.fa`, `.fna`, `.fasta`, `.gb`, `.gbk`, `.gbff`)
        outdir : str | Path | None, optional
            Temporary result directory. If None, tmp directory is auto created.
        threads : int | None, optional
            Threads parameter for lastal run
        cmd_opts : str | None, optional
            `lastal` additional command options
        """
        super().__init__()

        self._seqs = self._parse_input_gbk_and_fasta_seqs(seqs)
        self._outdir = None if outdir is None else Path(outdir)
        self._threads = self.max_threads if threads is None else threads
        self._cmd_opts = cmd_opts

    @classmethod
    def get_tool_name(cls) -> str:
        """Tool name"""
        return "Last"

    @classmethod
    def get_binary_names(cls) -> list[str]:
        """Binary names"""
        return ["lastdb", "lastal", "last-split", "maf-convert"]

    @classmethod
    def get_version(cls) -> str:
        """Tool version"""
        return cls._get_version(
            cmd="lastal --version",
            pattern=r"lastal (\S+)",
        )

    def run(self) -> list[AlignCoord]:
        """Run genome alignment"""
        # Last genome alignment parameters are selected based on this page
        # https://gitlab.com/mcfrith/last/-/blob/main/doc/last-cookbook.rst
        with TemporaryDirectory() as tmpdir:
            outdir = self._outdir if self._outdir else tmpdir
            outdir = Path(outdir)
            os.makedirs(outdir, exist_ok=True)
            genome_files: list[Path] = self._write_genome_files(self._seqs, outdir)

            logger.info(f"{'=' * 10} Start Last Alignment {'=' * 10}")
            align_coords = []
            for idx in range(len(genome_files) - 1):
                qfile, rfile = genome_files[idx], genome_files[idx + 1]
                qname, rname = qfile.stem, rfile.stem
                logger.info(f"{idx + 1:02d}. Last Alignment '{qname}' vs '{rname}'")
                # Make Last database
                lastdb = outdir / f"{rname}_lastdb"
                cmd = f"lastdb '{lastdb}' '{rfile}' -P {self._threads}"
                self.run_cmd(cmd)
                # Run Last alignment
                outfile_prefix = f"{idx + 1:02d}_{qname}_vs_{rname}"
                maf_outfile1 = outdir / f"{outfile_prefix}_many-to-one.maf"
                cmd = f"lastal '{lastdb}' '{qfile}' -P {self._threads} -D 1e9 --split-f=MAF+"  # noqa
                if self._cmd_opts:
                    cmd = f"{cmd} {self._cmd_opts}"
                self.run_cmd(cmd, maf_outfile1)
                # Convert many-to-one -> one-to-one
                maf_outfile2 = outdir / f"{outfile_prefix}_one-to-one.maf"
                cmd = f"last-split -r '{maf_outfile1}'"
                self.run_cmd(cmd, maf_outfile2)
                # Convert MAF -> BlastTab format
                blast_outfile = outdir / f"{outfile_prefix}.tsv"
                cmd = f"maf-convert -n blasttab '{maf_outfile2}'"
                self.run_cmd(cmd, blast_outfile)

                align_coords.extend(
                    AlignCoord.parse_blast_file(blast_outfile, qname, rname)
                )
            logger.info(f"{'=' * 10} Finish Last Alignment {'=' * 10}")

        return align_coords
