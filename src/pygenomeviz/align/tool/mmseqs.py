from __future__ import annotations

import csv
import logging
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Sequence

from pygenomeviz.align import AlignCoord
from pygenomeviz.align.tool import AlignToolBase
from pygenomeviz.parser import Genbank


class MMseqs(AlignToolBase):
    """MMseqs RBH Search Class"""

    def __init__(
        self,
        seqs: Sequence[str | Path | Genbank],
        *,
        outdir: str | Path | None = None,
        evalue: float = 1e-3,
        threads: int | None = None,
        cmd_opts: str | None = None,
        logger: logging.Logger | None = None,
        quiet: bool = True,
    ):
        """
        Parameters
        ----------
        seqs : Sequence[str | Path | Genbank]
            List of `genbank file` or `Genbank object`
        outdir : str | Path | None, optional
            Temporary result directory. If None, tmp directory is auto created.
        evalue : float, optional
            E-value parameter for MMseqs run
        threads : int | None, optional
            Threads parameter for MMseqs run
        cmd_opts : str | None, optional
            `mmseqs easy-rbh` additional command options
        logger : logging.Logger | None, optional
            Logger object. If None, logger instance newly created.
        quiet : bool, optional
            If True, don't display log message
        """
        super().__init__(logger, quiet)

        self._seqs = self._parse_input_gbk_seqs(seqs)
        self._outdir = None if outdir is None else Path(outdir)
        self._evalue = evalue
        self._threads = self.max_threads if threads is None else threads
        self._cmd_opts = cmd_opts

    @classmethod
    def get_tool_name(cls) -> str:
        """Tool name"""
        return "MMseqs"

    @classmethod
    def get_binary_names(cls) -> list[str]:
        """Binary names"""
        return ["mmseqs"]

    def run(self) -> list[AlignCoord]:
        """Run genome alignment"""
        with TemporaryDirectory() as tmpdir:
            outdir = self._outdir if self._outdir else tmpdir
            outdir = Path(outdir)
            os.makedirs(outdir, exist_ok=True)
            cds_files: list[Path] = self._write_cds_files(outdir)

            self._logger.info(f"{'='*10} Start MMseqs RBH Search {'='*10}")
            align_coords = []
            for idx in range(len(cds_files) - 1):
                qfile, rfile = cds_files[idx], cds_files[idx + 1]
                qname, rname = qfile.stem, rfile.stem
                rbh_file = outdir / f"{idx+1:02d}_{qname}_vs_{rname}.tsv"
                log_msg = f"{idx+1:02d}. MMseqs RBH Search '{qname}' vs '{rname}'"
                self._logger.info(log_msg)
                cmd = f"mmseqs easy-rbh '{qfile}' '{rfile}' '{rbh_file}' {outdir} --threads {self._threads} -e {self._evalue} -v 0"  # noqa: E501
                if self._cmd_opts:
                    cmd = f"{cmd} {self._cmd_opts}"
                self.run_cmd(cmd, self._logger)

                align_coords.extend(self._parse_coords_file(rbh_file, qname, rname))
            self._logger.info(f"{'='*10} Finish MMseqs RBH Search {'='*10}")

        return align_coords

    def _write_cds_files(self, outdir: str | Path) -> list[Path]:
        """Write CDS files"""
        outdir = Path(outdir)
        cds_files: list[Path] = []
        for gbk in self._seqs:
            cds_file = outdir / f"{gbk.name}.faa"
            self._logger.info(f"Convert Genbank object to CDS fasta file '{cds_file}'")
            gbk.write_cds_fasta(cds_file)
            cds_files.append(cds_file)
        return cds_files

    def _parse_coords_file(
        self,
        rbh_result_file: str | Path,
        query_id: str,
        ref_id: str,
    ) -> list[AlignCoord]:
        """Parse MMseqs RBH result

        Parameters
        ----------
        rbh_result_file : str | Path
            MMseqs RBH result file
        query_id : str
            Query ID
        ref_id : str
            Reference ID

        Returns
        -------
        align_coords : list[AlignCoord]
            Align Coords
        """
        align_coords = []
        dup_check_list = []
        with open(rbh_result_file) as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                # Query [e.g. `GENE000001_NC_XXXXXX.X|name|8205_8559_1|`]
                qname = row[0].split("|")[-3]
                qstart, qend, qstrand = [
                    int(d) for d in row[0].split("|")[-2].split("_")
                ]
                if qstrand == -1:
                    qstart, qend = qend, qstart
                # Reference
                rname = row[1].split("|")[-3]
                rstart, rend, rstrand = [
                    int(d) for d in row[1].split("|")[-2].split("_")
                ]
                if rstrand == -1:
                    rstart, rend = rend, rstart
                # Identity (%) [e.g. 0.95215 -> 95.21]
                identity = int(float(row[2]) * 10000) / 100
                evalue = float(row[10])

                # Check duplication
                qkey = f"{query_id} {qname} {qstart} {qend} {qstrand}"
                rkey = f"{ref_id} {rname} {rstart} {rend} {rstrand}"
                if rkey in dup_check_list:
                    continue
                if qkey in dup_check_list:
                    continue
                dup_check_list.extend([qkey, rkey])

                align_coord = AlignCoord(
                    query_id,
                    qname,
                    qstart,
                    qend,
                    ref_id,
                    rname,
                    rstart,
                    rend,
                    identity,
                    evalue,
                )
                align_coords.append(align_coord)

        return align_coords
