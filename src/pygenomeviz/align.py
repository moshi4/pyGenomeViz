from __future__ import annotations

import csv
import itertools
import multiprocessing as mp
import os
import platform
import shutil
import subprocess as sp
import sys
from dataclasses import astuple, dataclass
from pathlib import Path
from typing import ClassVar, List, Optional, Tuple, Union

from pygenomeviz import Genbank


class Align:
    """Run MUMmer Genome Alignment Class"""

    def __init__(
        self,
        genome_resources: Union[List[Union[str, Path]], List[Genbank]],
        outdir: Union[str, Path],
        seqtype: str = "nucleotide",
        maptype: str = "one-to-one",
        process_num: Optional[int] = None,
    ):
        """
        Parameters
        ----------
        genome_resources : Union[List[Union[str, Path]], List[Genbank]]
            Genome fasta files or Genbank objects
        outdir : Union[str, Path]
            Output directory
        seqtype : str, optional
            Sequence type (`nucleotide`|`protein`)
        maptype : str, optional
            Alignment map type (`one-to-one`|`many-to-many`)
        process_num : Optional[int], optional
            Use processor number (Default: `'Max Processor' - 1`)
        """
        self.genome_resources = genome_resources
        self.outdir = Path(outdir)
        self.seqtype = seqtype.lower()
        self.maptype = maptype.lower()
        self.process_num = self._max_process_num if process_num is None else process_num

        if not self.check_installation():
            err_msg = "ERROR: Genome alignment by MUMmer is not available "
            err_msg += "in this environment. Please check MUMmer installation."
            print(err_msg)
            sys.exit(1)

    @property
    def genome_num(self) -> int:
        """Input genome fasta file number"""
        return len(self.genome_resources)

    @property
    def _genome_fasta_files(self) -> List[Path]:
        """Genome fasta file list"""
        genome_fasta_files = []
        for gr in self.genome_resources:
            if isinstance(gr, Genbank):
                suffix = "_reverse.fna" if gr.reverse else ".fna"
                filename = f"{gr.name}_{gr.min_range}-{gr.max_range}{suffix}"
                genome_fasta_file = self.outdir / filename
                if not genome_fasta_file.exists():
                    gr.write_genome_fasta(genome_fasta_file)
                genome_fasta_files.append(genome_fasta_file)
            else:
                genome_fasta_files.append(Path(gr))
        return genome_fasta_files

    @property
    def _align_binary(self) -> str:
        """Genome alignment binary name (`nucleotide='nucmer'`|`protein='promer'`)"""
        if self.seqtype == "nucleotide":
            return "nucmer"
        elif self.seqtype == "protein":
            return "promer"
        else:
            raise ValueError(f"Invalid seqtype '{self.seqtype}'")

    @property
    def _map_option(self) -> str:
        """Genome alignment mapping option (`one-to-one='-1'`|`many-to-many='-m'`)"""
        if self.maptype == "one-to-one":
            return "-1"
        elif self.maptype == "many-to-many":
            return "-m"
        else:
            raise ValueError(f"Invalid maptype '{self.maptype}'")

    @property
    def _max_process_num(self) -> int:
        """Max process number"""
        cpu_num = os.cpu_count()
        return 1 if cpu_num is None or cpu_num == 1 else cpu_num - 1

    def run(self) -> List[AlignCoord]:
        """Run MUMmer genome alignment

        Returns
        -------
        align_coords : List[AlignCoord]
            Genome alignment coord list
        """
        # Prepare data for run MUMmer with multiprocessing
        mp_data_list: List[Tuple[Path, Path, int]] = []
        for idx in range(0, self.genome_num - 1):
            fa_file1 = self._genome_fasta_files[idx]
            fa_file2 = self._genome_fasta_files[idx + 1]
            mp_data_list.append((fa_file1, fa_file2, idx))

        # Run MUMmer with multiprocessing
        with mp.Pool(processes=self.process_num) as p:
            results = p.starmap(self._run_mummer, mp_data_list)

        return list(itertools.chain.from_iterable(results))

    def _run_mummer(self, fa_file1: Path, fa_file2: Path, idx: int) -> List[AlignCoord]:
        """Run MUMmer function for multiprocessing

        Parameters
        ----------
        fa_file1 : Path
            Genome fasta file 01
        fa_file2 : Path
            Genome fasta file 02
        idx : int
            Multiprocessing index

        Returns
        -------
        align_coords : List[AlignCoord]
            Align coord list
        """
        # Run genome alignment using nucmer or promer
        prefix = self.outdir / f"out{idx}"
        delta_file = prefix.with_suffix(".delta")
        cmd = f"{self._align_binary} {fa_file1} {fa_file2} --prefix={prefix}"
        _ = sp.run(cmd, shell=True, capture_output=True, text=True)

        # Run delta-filter to map 'one-to-one' or 'many-to-many' relation
        filter_delta_file = self.outdir / f"filter_out{idx}.delta"
        cmd = f"delta-filter {self._map_option} {delta_file} > {filter_delta_file}"
        _ = sp.run(cmd, shell=True, capture_output=True, text=True)

        # Run show-coords to extract alingment coords
        # -H: no header, -T: tab-delimited format, -r: sort by ref-id,
        # -k: knockout alignments that overlap in different frame (promer only)
        coords_file = self.outdir / f"coords{idx}.tsv"
        cmd = f"show-coords -H -T -r -k {filter_delta_file} > {coords_file}"
        _ = sp.run(cmd, shell=True, capture_output=True, text=True)

        align_coords = AlignCoord.parse(coords_file, self.seqtype)

        # Delete work files
        for work_file in (delta_file, filter_delta_file, coords_file):
            os.unlink(work_file)

        return align_coords

    @staticmethod
    def check_installation() -> bool:
        """Check MUMmer installation

        Returns
        -------
        result : bool
            Check result
        """
        # MacOS or Linux only
        if platform.system() not in ("Darwin", "Linux"):
            return False
        # Mummer binary installation check
        required_bins = ["nucmer", "promer", "delta-filter", "show-coords"]
        for required_bin in required_bins:
            if not shutil.which(required_bin):
                return False
        return True


@dataclass
class AlignCoord:
    """MUMmer Alignment Coordinates DataClass"""

    ref_start: int
    ref_end: int
    query_start: int
    query_end: int
    ref_length: int
    query_length: int
    identity: float
    ref_name: str
    query_name: str
    header_list: ClassVar[List[str]] = [
        "REF_START",
        "REF_END",
        "QUERY_START",
        "QUERY_END",
        "REF_LENGTH",
        "QUERY_LENGTH",
        "IDENTITY",
        "REF_NAME",
        "QUERY_NAME",
    ]

    @property
    def as_tsv_format(self) -> str:
        """TSV format text"""
        return "\t".join([str(v) for v in astuple(self)])

    @staticmethod
    def write(align_coords: List[AlignCoord], outfile: Union[str, Path]) -> None:
        """Write alignment results as tsv format file

        Parameters
        ----------
        align_coords : List[AlignCoord]
            Alignment coords
        outfile : Union[str, Path]
            Output file path
        """
        with open(outfile, "w") as f:
            header = "\t".join(AlignCoord.header_list)
            output = "\n".join([ac.as_tsv_format for ac in align_coords])
            f.write(header + "\n" + output)

    @staticmethod
    def parse(
        coords_tsv_file: Union[str, Path],
        seqtype: str,
    ) -> List[AlignCoord]:
        """Parse MUMmer(nucmer|promer) output coords result file

        Parameters
        ----------
        coords_tsv_file : Union[str, Path]
            MUMmer align coords file
        seqtype : str
            Sequence type (`nucleotide`|`protein`)

        Returns
        -------
        align_coords : List[AlignCoord]
            Align coord list
        """
        align_coords = []
        with open(coords_tsv_file) as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                # Check read file contents & extract required row values
                if seqtype == "nucleotide":
                    if len(row) != 9:
                        err_msg = f"Invalid nucmer coords file '{coords_tsv_file}'!!"
                        raise ValueError(err_msg)
                elif seqtype == "protein":
                    if len(row) != 13:
                        err_msg = f"Invalid promer coords file '{coords_tsv_file}'!!"
                        raise ValueError(err_msg)
                    row = row[0:7] + row[11:13]
                else:
                    raise ValueError(f"Invalid seqtype '{seqtype}'!!")

                # Convert to correct value type
                typed_row = []
                for idx, val in enumerate(row):
                    if 0 <= idx <= 5:
                        typed_row.append(int(val))
                    elif idx == 6:
                        typed_row.append(float(val))
                    else:
                        typed_row.append(str(val))

                align_coords.append(AlignCoord(*typed_row))

        return align_coords

    @staticmethod
    def filter(
        align_coords: List[AlignCoord],
        min_length: int = 0,
        min_identity: float = 0.0,
    ) -> List[AlignCoord]:
        """Filter align coord list with 'length' & 'identity'

        Parameters
        ----------
        align_coords : List[AlignCoord]
            Align coord list
        min_length : int, optional
            Min length filtering threshold
        min_identity : float, optional
            Min identity filtering threshold

        Returns
        -------
        filtered_align_coords : List[AlignCoord]
            Filtered align coord list
        """
        filtered_align_coords: List[AlignCoord] = []
        for ac in align_coords:
            rlen, qlen, ident = ac.ref_length, ac.query_length, ac.identity
            if (rlen >= min_length and qlen >= min_length) and ident >= min_identity:
                filtered_align_coords.append(AlignCoord(*astuple(ac)))
        return filtered_align_coords
