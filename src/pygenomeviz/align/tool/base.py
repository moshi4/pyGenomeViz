from __future__ import annotations

import logging
import os
import shlex
import shutil
import subprocess as sp
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Sequence

from pygenomeviz.align import AlignCoord
from pygenomeviz.logger import get_logger
from pygenomeviz.parser import Genbank


class AlignToolBase(ABC):
    """Alignment Tool Abstract Base Class"""

    def __init__(self, logger: logging.Logger | None, quiet: bool = False):
        """
        Parameters
        ----------
        logger : logging.Logger | None
            Logger object. If None, logger instance newly created.
        quiet : bool, optional
            If True, don't display log message.
        """
        self._logger = logger if logger else get_logger(__name__, quiet=quiet)
        self.check_installation(logger=self._logger)

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

    @abstractmethod
    def run(self) -> list[AlignCoord]:
        """Run genome alignment"""
        raise NotImplementedError

    @classmethod
    def check_installation(
        cls,
        exit_on_false: bool = True,
        logger: logging.Logger | None = None,
    ) -> bool:
        """Check required binaries installation

        Parameters
        ----------
        exit_on_false : bool, optional
            If True and check result is False, raise RuntimeError.
        logger : logging.Logger | None, optional
            Logger object. If None, logger instance newly created.

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
            logger = logger if logger else get_logger(__name__)
            logger.error(f"'{tool_name}' is not available in this environment.")
            logger.error(f"Please check '{tool_name}' installation.")
            raise RuntimeError(f"Failed to run '{tool_name}' aligner!!")

        return is_installed

    def run_cmd(
        self,
        cmd: str,
        logger: logging.Logger,
        stdout_file: str | Path | None = None,
    ) -> None:
        """Run command

        Parameters
        ----------
        cmd : str
            Command to run
        logger : logging.Logger
            Logger object
        stdout_file : str | Path | None, optional
            Write stdout result if file is set
        """
        logger.info(f"$ {cmd}")
        cmd_args = shlex.split(cmd)
        cmd_res = sp.run(cmd_args, capture_output=True, text=True)

        if cmd_res.returncode == 0:
            # Write stdout result if stdout_file is set
            if stdout_file:
                logger.info(f"> Save cmd stdout results to '{stdout_file}'")
                with open(stdout_file, "w") as f:
                    f.write(cmd_res.stdout)
        else:
            logger.error("Failed to run command below!!")
            logger.error(f"$ {cmd}")
            stdout_lines = cmd_res.stdout.splitlines()
            if len(stdout_lines) > 0:
                logger.error("STDOUT:")
                for line in stdout_lines:
                    logger.error(f"> {line}")
            stderr_lines = cmd_res.stderr.splitlines()
            if len(stderr_lines) > 0:
                logger.error("STDERR:")
                for line in stderr_lines:
                    logger.error(f"> {line}")
            raise RuntimeError(f"Failed to run '{self.get_tool_name()}' aligner!!")

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
        self, seqs: Sequence[str | Path | Genbank]
    ) -> Sequence[Path | Genbank]:
        """Parse input genbank and fasta sequences

        Parameters
        ----------
        seqs : Sequence[str | Path | Genbank]
            List of `fasta file` or `genbank file` or `Genbank object`

        Returns
        -------
        parse_seqs : Sequence[Path | Genbank]
            List of `fasta file` or `Genbank object`
        """
        # Check number of seqs
        if len(seqs) < 2:
            raise ValueError("Number of input seqs is less than 2.")
        # Parse genbank or fasta files
        parse_seqs: Sequence[Path | Genbank] = []
        for seq in seqs:
            if isinstance(seq, str):
                seq = Path(seq)
            if isinstance(seq, Path):
                gbk_suffixes = (".gb", ".gbk", ".gbff", ".gz")
                fasta_suffixes = (".fa", ".fna", ".fasta")
                if seq.suffix in gbk_suffixes:
                    parse_seqs.append(Genbank(seq))
                elif seq.suffix in fasta_suffixes:
                    parse_seqs.append(seq)
                else:
                    valid_suffixes = gbk_suffixes + fasta_suffixes
                    err_msg = f"'{seq}' is invalid file suffix ({valid_suffixes=})"
                    raise ValueError(err_msg)
            else:
                parse_seqs.append(seq)
        return parse_seqs
