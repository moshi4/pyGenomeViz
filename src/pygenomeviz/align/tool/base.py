from __future__ import annotations

import logging
import os
import re
import shlex
import shutil
import subprocess as sp
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Sequence

from pygenomeviz.align import AlignCoord
from pygenomeviz.const import UNKNOWN_VERSION
from pygenomeviz.parser import Fasta, Genbank

logger = logging.getLogger(__name__)


class AlignToolBase(ABC):
    """Alignment Tool Abstract Base Class"""

    def __init__(self):
        self.check_installation()

    @property
    def max_threads(self) -> int:
        """Max threads number"""
        cpu_num = os.cpu_count()
        return 1 if cpu_num is None or cpu_num == 1 else cpu_num

    @classmethod
    @abstractmethod
    def get_tool_name(cls) -> str:
        """Tool name"""
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def get_binary_names(cls) -> list[str]:
        """Binary names"""
        raise NotImplementedError

    @classmethod
    @abstractmethod
    def get_version(cls) -> str:
        """Tool version"""
        raise NotImplementedError

    @abstractmethod
    def run(self) -> list[AlignCoord]:
        """Run genome alignment"""
        raise NotImplementedError

    @classmethod
    def check_installation(
        cls,
        exit_on_false: bool = True,
    ) -> bool:
        """Check required binaries installation

        Parameters
        ----------
        exit_on_false : bool, optional
            If True and check result is False, raise RuntimeError.

        Returns
        -------
        result : bool
            Check result
        """
        is_installed = True
        for required_binary in cls.get_binary_names():
            if not shutil.which(required_binary):
                is_installed = False

        if not is_installed and exit_on_false:
            tool_name = cls.get_tool_name()
            logger.error(f"'{tool_name}' is not available in this environment.")
            logger.error(f"Please check '{tool_name}' installation.")
            raise RuntimeError(f"Failed to run '{tool_name}' aligner!!")

        return is_installed

    def run_cmd(
        self,
        cmd: str,
        stdout_file: str | Path | None = None,
    ) -> None:
        """Run command

        Parameters
        ----------
        cmd : str
            Command to run
        stdout_file : str | Path | None, optional
            Write stdout result if file is set
        """
        logger.info(f"$ {cmd}")
        cmd_args = shlex.split(cmd)
        try:
            cmd_res = sp.run(cmd_args, capture_output=True, text=True, check=True)
            # Write stdout result if stdout_file is set
            if stdout_file:
                logger.info(f"> Save cmd stdout results to '{stdout_file}'")
                with open(stdout_file, "w", encoding="utf-8") as f:
                    f.write(cmd_res.stdout)
        except sp.CalledProcessError as e:
            returncode, stdout, stderr = e.returncode, str(e.stdout), str(e.stderr)
            logger.error(f"Failed to run command below ({returncode=})")
            logger.error(f"$ {cmd}")
            stdout_lines = stdout.splitlines()
            if len(stdout_lines) > 0:
                logger.error("STDOUT:")
                for line in stdout_lines:
                    logger.error(f"> {line}")
            stderr_lines = stderr.splitlines()
            if len(stderr_lines) > 0:
                logger.error("STDERR:")
                for line in stderr_lines:
                    logger.error(f"> {line}")
            raise RuntimeError(f"Failed to run '{self.get_tool_name()}' aligner!!")

    @classmethod
    def _get_version(cls, cmd: str, pattern: str) -> str:
        """Get tool version by cmd & regex pattern

        Parameters
        ----------
        cmd : str
            Command to get version info
        pattern : str
            Regex pattern for search version

        Returns
        -------
        version : str
            Tool version (e.g. `v1.2.3`)
        """
        try:
            cmd_args = shlex.split(cmd)
            cmd_res = sp.run(cmd_args, capture_output=True, text=True)
            output = cmd_res.stderr if cmd_res.stdout == "" else cmd_res.stdout
            version = re.findall(pattern, output, re.MULTILINE)[0]
            return version
        except Exception:
            return UNKNOWN_VERSION

    def _parse_input_gbk_seqs(
        self, seqs: Sequence[str | Path | Genbank]
    ) -> list[Genbank]:
        """Parse input genbank sequences

        Parameters
        ----------
        seqs : Sequence[str | Path | Genbank]
            List of `genbank file` or `Genbank object`

        Returns
        -------
        parse_seqs : list[Genbank]
            List of `Genbank object`
        """
        # Check number of seqs
        if len(seqs) < 2:
            raise ValueError("Number of input seqs is less than 2.")
        # Parse genbank
        parse_seqs: list[Genbank] = []
        for seq in seqs:
            if isinstance(seq, (str, Path)):
                parse_seqs.append(Genbank(seq))
            else:
                parse_seqs.append(seq)
        return parse_seqs

    def _parse_input_gbk_and_fasta_seqs(
        self, seqs: Sequence[str | Path | Fasta | Genbank]
    ) -> Sequence[Fasta | Genbank]:
        """Parse input genbank and fasta sequences

        Parameters
        ----------
        seqs : Sequence[str | Path | Fasta | Genbank]
            List of fasta or genbank

        Returns
        -------
        parse_seqs : Sequence[Fasta | Genbank]
            List of fasta or genbank
        """
        # Check number of seqs
        if len(seqs) < 2:
            raise ValueError("Number of input seqs is less than 2.")
        # Parse genbank or fasta files
        parse_seqs: Sequence[Fasta | Genbank] = []
        gbk_suffixes = (".gb", ".gbk", ".gbff", ".gb.gz", ".gbk.gz", ".gbff.gz")
        fasta_suffixes = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")
        for seq in seqs:
            if isinstance(seq, str):
                seq = Path(seq)
            if isinstance(seq, Path):
                suffix = seq.suffix
                if suffix == ".gz" and len(seq.suffixes) >= 2:
                    suffix = "".join(seq.suffixes[-2])
                if suffix in gbk_suffixes:
                    parse_seqs.append(Genbank(seq))
                elif suffix in fasta_suffixes:
                    parse_seqs.append(Fasta(seq))
                else:
                    valid_suffixes = gbk_suffixes + fasta_suffixes
                    raise ValueError(f"'{seq}' is invalid file suffix ({valid_suffixes=})")  # fmt: skip  # noqa: E501
            else:
                parse_seqs.append(seq)
        return parse_seqs

    def _write_genome_files(
        self,
        seqs: Sequence[Fasta | Genbank],
        outdir: str | Path,
    ) -> list[Path]:
        """Write genome fasta files to output directory

        Parameters
        ----------
        seqs : Sequence[Fasta | Genbank]
            List of fasta or genbank
        outdir : str | Path
            Target output directory

        Returns
        -------
        genome_files : list[Path]
            Genome fasta files
        """
        genome_files: list[Path] = []
        for seq in seqs:
            genome_file = Path(outdir) / f"{seq.name}.fna"
            cls_name = seq.__class__.__name__
            logger.info(f"Convert {cls_name} object to genome fasta file '{genome_file}'")  # fmt: skip  # noqa: E501
            seq.write_genome_fasta(genome_file)
            genome_files.append(genome_file)
        return genome_files
