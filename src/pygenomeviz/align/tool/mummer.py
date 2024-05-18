from __future__ import annotations

import logging
import os
import shutil
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Literal, Sequence, get_args

from pygenomeviz.align import AlignCoord
from pygenomeviz.align.tool import AlignToolBase
from pygenomeviz.parser import Genbank
from pygenomeviz.typing import SeqType


class MUMmer(AlignToolBase):
    """MUMmer Alignment Class"""

    def __init__(
        self,
        seqs: Sequence[str | Path | Genbank],
        *,
        outdir: str | Path | None = None,
        seqtype: SeqType = "nucleotide",
        maptype: Literal["one-to-one", "many-to-many"] = "one-to-one",
        threads: int | None = None,
        cmd_opts: str | None = None,
        logger: logging.Logger | None = None,
        quiet: bool = True,
    ):
        """
        Parameters
        ----------
        seqs : Sequence[str | Path | Genbank]
            List of `fasta file` or `genbank file` or `Genbank object`
            (file suffix must be `.fa`, `.fna`, `.fasta`, `.gb`, `.gbk`, `.gbff`)
        outdir : str | Path | None, optional
            Temporary result directory. If None, tmp directory is auto created.
        seqtype : SeqType, optional
            `nucleotide`(blastn) or `protein`(tblastx)
        maptype : str, optional
            `one-to-one` or `many-to-many`
        evalue : float, optional
            E-value parameter for blast run
        threads : int | None, optional
            Threads parameter for blast run
        cmd_opts : str | None, optional
            `nucmer` or `promer` additional command options
        logger : logging.Logger | None, optional
            Logger object. If None, logger instance newly created.
        quiet : bool, optional
            If True, don't display log message
        """
        super().__init__(logger, quiet)

        valid_seqtype = get_args(SeqType)
        if seqtype not in valid_seqtype:
            raise ValueError(f"{seqtype=} is invalid ({valid_seqtype=})")

        self._seqs = self._parse_input_gbk_and_fasta_seqs(seqs)
        self._outdir = None if outdir is None else Path(outdir)
        self._seqtype: SeqType = seqtype
        self._maptype = maptype
        self._threads = self.max_threads if threads is None else threads
        self._cmd_opts = cmd_opts

    @classmethod
    def get_tool_name(cls) -> str:
        """Tool name"""
        return "MUMmer"

    @classmethod
    def get_binary_names(cls) -> list[str]:
        """Binary names"""
        return ["nucmer", "promer", "delta-filter", "show-coords"]

    def run(self) -> list[AlignCoord]:
        """Run genome alignment"""
        with TemporaryDirectory() as tmpdir:
            # Genome files
            outdir = self._outdir if self._outdir else tmpdir
            outdir = Path(outdir)
            os.makedirs(outdir, exist_ok=True)
            genome_files: list[Path] = self._write_genome_files(outdir)

            # Run MUMmer(nucmer/promer)
            self._logger.info(f"{'='*10} Start MUMmer Alignment {'='*10}")
            align_coords: list[AlignCoord] = []
            for idx in range(len(genome_files) - 1):
                qfile, rfile = genome_files[idx], genome_files[idx + 1]
                qname, rname = qfile.stem, rfile.stem
                log_msg = f"{idx+1:02d}. MUMmer Alignment '{qname}' vs '{rname}'"
                self._logger.info(log_msg)

                # Run genome alignment using nucmer or promer
                prefix = outdir / f"{idx+1:02d}_{qname}_vs_{rname}"
                seqtype2align_tool = dict(nucleotide="nucmer", protein="promer")
                align_tool = seqtype2align_tool[self._seqtype]
                cmd = f"{align_tool} --mum '{rfile}' '{qfile}' --prefix={prefix}"
                if self._cmd_opts:
                    cmd = f"{cmd} {self._cmd_opts}"
                self.run_cmd(cmd, self._logger)

                # Run delta-filter to map 'one-to-one' or 'many-to-many' relation
                delta_file = Path(str(prefix) + ".delta")
                filter_delta_file = Path(str(prefix) + "_filter.delta")
                maptype2mapopt = {"one-to-one": "-1", "many-to-many": "-m"}
                mapopt = maptype2mapopt[self._maptype]
                cmd = f"delta-filter {mapopt} '{delta_file}'"
                self.run_cmd(cmd, self._logger, filter_delta_file)

                # Run show-coords to extract alignment coords
                # -H: no header, -T: tab-delimited format, -r: sort by ref-id,
                # -k: knockout alignments that overlap in different frame (promer only)
                coords_file = Path(str(prefix) + "_coords.tsv")
                cmd = f"show-coords -H -T -r -k '{filter_delta_file}'"
                self.run_cmd(cmd, self._logger, coords_file)

                align_coords.extend(
                    AlignCoord.parse_mummer_file(
                        coords_file, qname, rname, self._seqtype
                    )
                )
            self._logger.info(f"{'='*10} Finish MUMmer Alignment {'='*10}")

        return align_coords

    def _write_genome_files(self, outdir: str | Path) -> list[Path]:
        """Write (or copy) genome fasta files to output directory

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
            if isinstance(seq, Genbank):
                genome_file = Path(outdir) / f"{seq.name}.fna"
                log_msg = f"Convert Genbank object to genome fasta file '{genome_file}'"
                self._logger.info(log_msg)
                seq.write_genome_fasta(genome_file)
                genome_files.append(genome_file)
            else:
                genome_file = Path(outdir) / seq.name
                self._logger.info(f"Copy fasta file to '{genome_file}'")
                shutil.copy(seq, genome_file)
                genome_files.append(genome_file)
        return genome_files
