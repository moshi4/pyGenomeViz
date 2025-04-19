from __future__ import annotations

import logging
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Sequence, get_args

from pygenomeviz.align import AlignCoord
from pygenomeviz.align.tool import AlignToolBase
from pygenomeviz.parser import Fasta, Genbank
from pygenomeviz.typing import SeqType

logger = logging.getLogger(__name__)


class Blast(AlignToolBase):
    """Blast Alignment Class"""

    def __init__(
        self,
        seqs: Sequence[str | Path | Fasta | Genbank],
        *,
        outdir: str | Path | None = None,
        seqtype: SeqType = "nucleotide",
        evalue: float = 1e-3,
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
        seqtype : SeqType, optional
            `nucleotide`(blastn) or `protein`(tblastx)
        evalue : float, optional
            E-value parameter for blast run
        threads : int | None, optional
            Threads parameter for blast run
        cmd_opts : str | None, optional
            `blastn` or `tblastx` additional command options
        """
        super().__init__()

        valid_seqtype = get_args(SeqType)
        if seqtype not in valid_seqtype:
            raise ValueError(f"{seqtype=} is invalid ({valid_seqtype=})")

        self._seqs = self._parse_input_gbk_and_fasta_seqs(seqs)
        self._outdir = None if outdir is None else Path(outdir)
        self._seqtype = seqtype
        self._evalue = evalue
        self._threads = self.max_threads if threads is None else threads
        self._cmd_opts = cmd_opts

    @classmethod
    def get_tool_name(cls) -> str:
        """Tool name"""
        return "BLAST"

    @classmethod
    def get_binary_names(cls) -> list[str]:
        """Binary names"""
        return ["makeblastdb", "blastn", "tblastx"]

    @classmethod
    def get_version(cls) -> str:
        """Tool version"""
        return cls._get_version(
            cmd="makeblastdb -version",
            pattern=r"makeblastdb: (\d+.\d+.\d+)",
        )

    def run(self) -> list[AlignCoord]:
        """Run genome alignment"""
        with TemporaryDirectory() as tmpdir:
            outdir = self._outdir if self._outdir else tmpdir
            outdir = Path(outdir)
            os.makedirs(outdir, exist_ok=True)
            genome_files: list[Path] = self._write_genome_files(self._seqs, outdir)

            logger.info(f"{'=' * 10} Start Blast Search {'=' * 10}")
            align_coords = []
            for idx in range(len(genome_files) - 1):
                qfile, rfile = genome_files[idx], genome_files[idx + 1]
                qname, rname = qfile.stem, rfile.stem
                logger.info(f"{idx + 1:02d}. Blast Search '{qname}' vs '{rname}'")
                # Make blast database
                blastdb = outdir / f"{rname}_blastdb"
                cmd = f"makeblastdb -in '{rfile}' -dbtype nucl -out '{blastdb}'"
                self.run_cmd(cmd)
                # Blast search ('blastn' or 'tblastx')
                seqtype2blast_tool = dict(nucleotide="blastn", protein="tblastx")
                blast_tool = seqtype2blast_tool[self._seqtype]
                blast_outfile = outdir / f"{idx + 1:02d}_{qname}_vs_{rname}.tsv"
                cmd = f"{blast_tool} -query '{qfile}' -db '{blastdb}' -out '{blast_outfile}' -outfmt 6 -evalue {self._evalue} -num_threads {self._threads}"  # noqa: E501
                if self._cmd_opts:
                    cmd = f"{cmd} {self._cmd_opts}"
                self.run_cmd(cmd)

                align_coords.extend(
                    AlignCoord.parse_blast_file(blast_outfile, qname, rname)
                )
            logger.info(f"{'=' * 10} Finish Blast Search {'=' * 10}")

        return align_coords
