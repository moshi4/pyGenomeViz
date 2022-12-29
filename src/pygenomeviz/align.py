from __future__ import annotations

import csv
import itertools
import multiprocessing as mp
import os
import shutil
import subprocess as sp
import sys
import tempfile
from abc import ABCMeta, abstractmethod
from dataclasses import astuple, dataclass
from pathlib import Path
from typing import ClassVar, Literal

from pygenomeviz import Genbank


class AlignToolBase(metaclass=ABCMeta):
    """Alignment Tool Abstract Base Class"""

    # Program name
    NAME: str = ""
    # Required binary names to run
    BINARIES: list[str] = []

    @abstractmethod
    def run(self) -> list[AlignCoord]:
        """Run genome alignment

        Returns
        -------
        align_coords : list[AlignCoord]
            Genome alignment coord list
        """
        raise NotImplementedError

    @property
    def max_process_num(self) -> int:
        """Max process number"""
        cpu_num = os.cpu_count()
        return 1 if cpu_num is None or cpu_num == 1 else cpu_num - 1

    @classmethod
    def check_installation(cls, exit_on_false: bool = True) -> bool:
        """Check required binaries installation

        Parameters
        ----------
        exit_on_error : bool
            If True and check result is False, system exit

        Returns
        -------
        result : bool
            Check result
        """
        is_installed = True
        for required_binary in cls.BINARIES:
            if not shutil.which(required_binary):
                is_installed = False

        if not is_installed and exit_on_false:
            err_msg = f"ERROR: Genome alignment by {cls.NAME} is not available "
            err_msg += f"in this environment. Please check {cls.NAME} installation."
            print(err_msg)
            sys.exit(1)

        return is_installed


class MUMmer(AlignToolBase):
    """MUMmer Alignment Class"""

    NAME = "MUMmer"
    BINARIES = ["nucmer", "promer", "delta-filter", "show-coords"]

    def __init__(
        self,
        gbk_resources: list[Genbank],
        outdir: str | Path,
        seqtype: Literal["protein", "nucleotide"] = "protein",
        maptype: Literal["many-to-many", "one-to-one"] = "many-to-many",
        process_num: int | None = None,
    ):
        """
        Parameters
        ----------
        gbk_resources : list[Genbank]
            Genbank objects
        outdir : str | Path
            Output directory
        seqtype : str, optional
            Sequence type (`nucleotide`|`protein`)
        maptype : str, optional
            Alignment map type (`one-to-one`|`many-to-many`)
        process_num : int | None, optional
            Use processor number (Default: `'Max Processor' - 1`)
        """
        self.check_installation()

        self.gbk_resources = gbk_resources
        self.outdir = Path(outdir)
        self.seqtype = seqtype
        self.maptype = maptype
        self.process_num = self.max_process_num if process_num is None else process_num

    @property
    def genome_num(self) -> int:
        """Input genome fasta file number"""
        return len(self.gbk_resources)

    @property
    def _genome_fasta_files(self) -> list[Path]:
        """Genome fasta file list"""
        genome_fasta_files = []
        for gr in self.gbk_resources:
            suffix = "_reverse.fna" if gr.reverse else ".fna"
            filename = f"{gr.name}_{gr.min_range}-{gr.max_range}{suffix}"
            genome_fasta_file = self.outdir / filename
            if not genome_fasta_file.exists():
                gr.write_genome_fasta(genome_fasta_file)
            genome_fasta_files.append(genome_fasta_file)
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

    def run(self) -> list[AlignCoord]:
        """Run genome alignment

        Returns
        -------
        align_coords : list[AlignCoord]
            Genome alignment coord list
        """
        # Prepare data for run MUMmer with multiprocessing
        mp_data_list: list[tuple[Path, Path, int]] = []
        for idx in range(0, self.genome_num - 1):
            fa_file1 = self._genome_fasta_files[idx]
            fa_file2 = self._genome_fasta_files[idx + 1]
            mp_data_list.append((fa_file1, fa_file2, idx))

        # Run MUMmer with multiprocessing
        with mp.Pool(processes=self.process_num) as p:
            results = p.starmap(self._run_mummer, mp_data_list)

        align_coords = list(itertools.chain.from_iterable(results))

        # Add min_range to start, end coordinates
        gbk_name2min_range = {gbk.name: gbk.min_range for gbk in self.gbk_resources}
        for ac in align_coords:
            ac.ref_start = ac.ref_start + gbk_name2min_range[ac.ref_name]
            ac.ref_end = ac.ref_end + gbk_name2min_range[ac.ref_name]
            ac.query_start = ac.query_start + gbk_name2min_range[ac.query_name]
            ac.query_end = ac.query_end + gbk_name2min_range[ac.query_name]

        return align_coords

    def _run_mummer(self, fa_file1: Path, fa_file2: Path, idx: int) -> list[AlignCoord]:
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
        align_coords : list[AlignCoord]
            Align coord list
        """
        # Run genome alignment using nucmer or promer
        prefix = self.outdir / f"out{idx}"
        delta_file = prefix.with_suffix(".delta")
        cmd = f"{self._align_binary} --mum {fa_file1} {fa_file2} --prefix={prefix}"
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

        align_coords = self.parse_coords_file(coords_file, self.seqtype)

        # Delete work files
        for work_file in (delta_file, filter_delta_file, coords_file):
            os.unlink(work_file)

        return align_coords

    @staticmethod
    def parse_coords_file(
        coords_tsv_file: str | Path,
        seqtype,
    ) -> list[AlignCoord]:
        """Parse MUMmer(nucmer|promer) output coords result file

        Parameters
        ----------
        coords_tsv_file : str | Path
            MUMmer align coords file

        Returns
        -------
        align_coords : list[AlignCoord]
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


class MMseqs(AlignToolBase):
    """MMseqs Alignment Class"""

    NAME = "MMseqs"
    BINARIES = ["mmseqs"]

    def __init__(
        self,
        gbk_resources: list[str | Path] | list[Genbank],
        outdir: str | Path,
        identity: float = 0,
        evalue: float = 1e-3,
        process_num: int | None = None,
    ):
        """
        Parameters
        ----------
        gbk_resources : list[str | Path] | list[Genbank]
            Genome sequence genbank files or Genbank objects
        outdir : str | Path
            Output directory
        identity : float, optional
            Identity threshold
        evalue : float, optional
            E-value threshold
        process_num : int | None, optional
            Use processor number (Default: `'Max Processor' - 1`)
        """
        self.check_installation()

        self.gbk_list: list[Genbank] = []
        for gr in gbk_resources:
            if isinstance(gr, Genbank):
                self.gbk_list.append(gr)
            else:
                self.gbk_list.append(Genbank(gr))
        self.outdir = Path(outdir)
        self.identity = identity
        self.evalue = evalue
        self.process_num = self.max_process_num if process_num is None else process_num

        os.makedirs(self.outdir, exist_ok=True)

    def run(self) -> list[AlignCoord]:
        """Run genome alignment

        Returns
        -------
        align_coords : list[AlignCoord]
            Genome alignment coord list
        """
        # Make CDS fasta files from Genbank files
        cds_fasta_files: list[Path] = []
        for gbk in self.gbk_list:
            cds_fasta_file = self.outdir / (gbk.name + ".faa")
            gbk.write_cds_fasta(cds_fasta_file)
            cds_fasta_files.append(cds_fasta_file)

        align_coords = []
        for idx in range(0, len(self.gbk_list) - 1):
            fa_file1, fa_file2 = cds_fasta_files[idx], cds_fasta_files[idx + 1]

            with tempfile.TemporaryDirectory() as tmpdir:
                name1 = fa_file1.with_suffix("").name
                name2 = fa_file2.with_suffix("").name
                rbh_result_file = self.outdir / f"{idx+1:02d}_{name1}-{name2}_rbh.tsv"
                cmd = f"mmseqs easy-rbh {fa_file1} {fa_file2} {rbh_result_file} "
                cmd += f"{tmpdir} --threads {self.process_num} -e {self.evalue} -v 0"
                print(f"# {idx+1:02d}: {name1}-{name2} RBH search\n$ {cmd}\n")
                sp.run(cmd, shell=True)

                align_coords.extend(
                    self.parse_rbh_result(rbh_result_file, name1, name2)
                )

        return align_coords

    def parse_rbh_result(
        self,
        rbh_result_file: str | Path,
        ref_name: str,
        query_name: str,
    ) -> list[AlignCoord]:
        """Parse MMseqs RBH result

        Parameters
        ----------
        rbh_result_file : str | Path
            MMseqs RBH result file
        ref_name : str
            Reference name
        query_name : str
            Query name

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
                # Reference
                ref_start, ref_end, ref_strand = [
                    int(d) for d in row[0].split("|")[1].split("_")
                ]
                ref_length = ref_end - ref_start + 1
                if ref_strand == -1:
                    ref_start, ref_end = ref_end, ref_start
                # Query
                query_start, query_end, query_strand = [
                    int(d) for d in row[1].split("|")[1].split("_")
                ]
                query_length = query_end - query_start + 1
                if query_strand == -1:
                    query_start, query_end = query_end, query_start
                # Identity (%) [e.g. 0.95215 -> 95.21]
                identity = int(float(row[2]) * 10000) / 100

                # Check duplication
                ref_key = f"{ref_name} {ref_start} {ref_end} {ref_strand}"
                query_key = f"{query_name} {query_start} {query_end} {query_strand}"
                if ref_key in dup_check_list:
                    continue
                if query_key in dup_check_list:
                    continue
                dup_check_list.extend([ref_key, query_key])

                align_coord = AlignCoord(
                    ref_start,
                    ref_end,
                    query_start,
                    query_end,
                    ref_length,
                    query_length,
                    identity,
                    ref_name,
                    query_name,
                )
                align_coords.append(align_coord)

        return align_coords


class ProgressiveMauve(AlignToolBase):
    """progresiveMauve Alignment Class"""

    NAME = "progressiveMauve"
    BINARIES = ["progressiveMauve"]

    def __init__(
        self,
        seq_files: list[str | Path],
        outdir: str | Path,
        refid: int = 0,
    ):
        """
        Parameters
        ----------
        seq_files : list[str | Path]
            Genome sequence files (Genbank or Fasta format)
        outdir : str | Path
            Output directory
        refid : int, optional
            Reference genome index
        """
        self.check_installation()

        self.seq_files = [Path(f) for f in seq_files]
        self.outdir = Path(outdir)
        self.refid = refid

        self.xmfa_file = self.outdir / "mauve.xmfa"
        self.bbone_file = self.outdir / "mauve_bbone.tsv"

        self.seq_outdir = self.outdir / "seqfiles"
        os.makedirs(self.seq_outdir, exist_ok=True)

    @property
    def filenames(self) -> list[str]:
        """File names"""
        return [f.with_suffix("").name for f in self.seq_files]

    def run(self) -> list[AlignCoord]:
        """Run genome alignment

        Returns
        -------
        align_coords : list[AlignCoord]
            Genome alignment coord list
        """
        # Copy seqfiles to output directory
        seq_copy_files = []
        for seq_file in self.seq_files:
            seq_copy_file = self.seq_outdir / Path(seq_file).name
            # progressiveMauve cannot recognize *.gbff as genbank format
            if str(seq_copy_file).endswith(".gbff"):
                seq_copy_file = seq_copy_file.with_suffix(".gbk")
            shutil.copy(seq_file, seq_copy_file)
            seq_copy_files.append(str(seq_copy_file))

        # Run progressiveMauve
        cmd = f"progressiveMauve --output={self.xmfa_file} "
        cmd += f"--backbone-output={self.bbone_file} {' '.join(seq_copy_files)}"
        sp.run(cmd, shell=True)

        return self.parse_pmauve_file(self.bbone_file)

    def parse_pmauve_file(self, bbone_file: str | Path) -> list[AlignCoord]:
        """Parse progressiveMauve bbone file

        Parameters
        ----------
        bbone_file : str | Path
            progressiveMauve bbone format file

        Returns
        -------
        align_coords : list[AlignCoord]
            Genome alignment coord list
        """
        with open(bbone_file) as f:
            reader = csv.reader(f, delimiter="\t")
            header_row = next(reader)
            genome_num = int(len(header_row) / 2)
            rows = []
            for row in reader:
                row = [int(col) for col in row]
                ref_idx = self.refid * 2
                # Always set reference seq coordinates to positive value
                if row[ref_idx] < 0:
                    row = [col * -1 for col in row]
                # Ignore no commonly conserved regions in all genomes
                if row.count(0) >= 2:
                    continue
                # Ignore too short conserved regions (< 20bp)
                if abs(row[ref_idx] - row[ref_idx + 1]) < 20:
                    continue
                rows.append(row)
            # Sort by reference seq coordinates
            rows = sorted(rows, key=lambda row: row[ref_idx])

        align_coords = []
        for row in rows:
            for i in range(genome_num - 1):
                idx = i * 2
                # Reference start-end
                rstart, rend = row[idx], row[idx + 1]
                if rstart < 0 and rend < 0:
                    rstart, rend = abs(rend), abs(rstart)
                # Query start-end
                qstart, qend = row[idx + 2], row[idx + 3]
                if qstart < 0 and qend < 0:
                    qstart, qend = abs(qend), abs(qstart)
                rlength, qlength = abs(rend - rstart) + 1, abs(qend - qstart) + 1
                rname, qname = self.filenames[i], self.filenames[i + 1]
                align_coord = AlignCoord(
                    rstart, rend, qstart, qend, rlength, qlength, 0, rname, qname
                )
                align_coords.append(align_coord)
        return align_coords


@dataclass
class AlignCoord:
    """Alignment Coordinates DataClass"""

    ref_start: int
    ref_end: int
    query_start: int
    query_end: int
    ref_length: int
    query_length: int
    identity: float
    ref_name: str
    query_name: str
    header_list: ClassVar[list[str]] = [
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
    def ref_strand(self) -> int:
        """Reference strand"""
        return 1 if self.ref_end > self.ref_start else -1

    @property
    def ref_link(self) -> tuple[str, int, int]:
        """Reference (name, start, end) link"""
        return (self.ref_name, self.ref_start, self.ref_end)

    @property
    def ref_block(self) -> tuple[int, int, int]:
        """Reference (start, end, strand) block"""
        if self.ref_start < self.ref_end:
            return (self.ref_start, self.ref_end, self.ref_strand)
        else:
            return (self.ref_end, self.ref_start, self.ref_strand)

    @property
    def query_strand(self) -> int:
        """Query strand"""
        return 1 if self.query_end > self.query_start else -1

    @property
    def query_link(self) -> tuple[str, int, int]:
        """Query (name, start, end) link"""
        return (self.query_name, self.query_start, self.query_end)

    @property
    def query_block(self) -> tuple[int, int, int]:
        """Query (start, end, strand) block"""
        if self.query_start < self.query_end:
            return (self.query_start, self.query_end, self.query_strand)
        else:
            return (self.query_end, self.query_start, self.query_strand)

    @property
    def is_inverted(self) -> bool:
        """Check inverted or not"""
        return self.ref_strand * self.query_strand < 0

    @property
    def as_tsv_format(self) -> str:
        """TSV format text"""
        return "\t".join([str(v) for v in astuple(self)])

    @staticmethod
    def write(align_coords: list[AlignCoord], outfile: str | Path) -> None:
        """Write alignment coords as tsv format file

        Parameters
        ----------
        align_coords : list[AlignCoord]
            Alignment coords
        outfile : str | Path
            Output file path
        """
        with open(outfile, "w") as f:
            header = "\t".join(AlignCoord.header_list)
            output = "\n".join([ac.as_tsv_format for ac in align_coords])
            f.write(header + "\n" + output)

    @staticmethod
    def read(align_coords_file: str | Path) -> list[AlignCoord]:
        """Read alignment coords tsv format file

        Parameters
        ----------
        align_coords_file : str | Path
            Alignment coords tsv file

        Returns
        -------
        align_coords : list[AlignCoord]
            Alignment coords
        """
        align_coords = []
        with open(align_coords_file) as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader)
            for row in reader:
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
        align_coords: list[AlignCoord],
        min_length: int = 0,
        min_identity: float = 0.0,
    ) -> list[AlignCoord]:
        """Filter align coord list with 'length' & 'identity'

        Parameters
        ----------
        align_coords : list[AlignCoord]
            Align coord list
        min_length : int, optional
            Min length filtering threshold
        min_identity : float, optional
            Min identity filtering threshold

        Returns
        -------
        filtered_align_coords : list[AlignCoord]
            Filtered align coord list
        """
        filtered_align_coords: list[AlignCoord] = []
        for ac in align_coords:
            rlen, qlen, ident = ac.ref_length, ac.query_length, ac.identity
            if (rlen >= min_length and qlen >= min_length) and ident >= min_identity:
                filtered_align_coords.append(AlignCoord(*astuple(ac)))
        return filtered_align_coords
