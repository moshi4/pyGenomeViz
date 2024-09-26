from __future__ import annotations

import bz2
import gzip
import zipfile
from io import TextIOWrapper
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class Fasta:
    """Fasta Parser Class"""

    def __init__(
        self,
        fasta: str | Path,
        *,
        name: str | None = None,
    ):
        """
        Parameters
        ----------
        fasta : str | Path
            Fasta file
        name : str | None, optional
            name (If None, file name is set)
        """
        fasta = Path(fasta)
        self._records = self._parse_fasta_file(fasta)

        # Set fasta name
        if name is not None:
            self._name = name
        else:
            if fasta.suffix in (".gz", ".bz2", ".zip"):
                self._name = fasta.with_suffix("").with_suffix("").name
            else:
                self._name = fasta.with_suffix("").name

        if len(self.records) == 0:
            raise ValueError(f"Failed to parse '{fasta}' as fasta file.")
        if len(self.records) != len(self.get_seqid2seq()):
            raise ValueError("Duplicate IDs are contained in fasta file.")

    ############################################################
    # Property
    ############################################################

    @property
    def name(self) -> str:
        """Name"""
        return self._name

    @property
    def records(self) -> list[SeqRecord]:
        """Fasta records"""
        return self._records

    @property
    def genome_seq(self) -> str:
        """Genome sequence (only first record)"""
        return str(self.records[0].seq)

    @property
    def genome_length(self) -> int:
        """Genome length (only first record)"""
        return len(self.genome_seq)

    @property
    def full_genome_seq(self) -> str:
        """Full genome sequence (concatenate all records)"""
        return "".join(str(r.seq) for r in self.records)

    @property
    def full_genome_length(self) -> int:
        """Full genome length (concatenate all records)"""
        return len(self.full_genome_seq)

    ############################################################
    # Public Method
    ############################################################

    def get_seqid2seq(self) -> dict[str, str]:
        """Get seqid & complete/contig/scaffold genome sequence dict

        Returns
        -------
        seqid2seq : dict[str, str]
            seqid & genome sequence dict
        """
        return {str(rec.id): str(rec.seq) for rec in self.records}

    def get_seqid2size(self) -> dict[str, int]:
        """Get seqid & complete/contig/scaffold genome size dict

        Returns
        -------
        seqid2size : dict[str, int]
            seqid & genome size dict
        """
        return {seqid: len(seq) for seqid, seq in self.get_seqid2seq().items()}

    def get_seqid2record(self) -> dict[str, SeqRecord]:
        """Get seqid & complete/contig/scaffold genome record dict

        Returns
        -------
        seqid2record : dict[str, SeqRecord]
            seqi & genome record dict
        """
        return {str(rec.id): rec for rec in self.records}

    def write_genome_fasta(self, outfile: str | Path) -> None:
        """Write genome fasta file

        Parameters
        ----------
        outfile : str | Path
            Output genome fasta file
        """
        with open(outfile, "w", encoding="utf-8") as f:
            for seqid, seq in self.get_seqid2seq().items():
                f.write(f">{seqid}\n{seq}\n")

    ############################################################
    # Private Method
    ############################################################

    def _parse_fasta_file(self, fasta_file: str | Path) -> list[SeqRecord]:
        """Parse fasta file

        Parameters
        ----------
        fasta_file : str | Path
            Fasta file

        Returns
        -------
        seq_records : list[SeqRecord]
            SeqRecord list
        """
        if Path(fasta_file).suffix == ".gz":
            with gzip.open(fasta_file, mode="rt", encoding="utf-8") as f:
                return list(SeqIO.parse(f, "fasta"))
        elif Path(fasta_file).suffix == ".bz2":
            with bz2.open(fasta_file, mode="rt", encoding="utf-8") as f:
                return list(SeqIO.parse(f, "fasta"))
        elif Path(fasta_file).suffix == ".zip":
            with zipfile.ZipFile(fasta_file) as zip:
                with zip.open(zip.namelist()[0]) as f:
                    io = TextIOWrapper(f, encoding="utf-8")
                    return list(SeqIO.parse(io, "fasta"))
        else:
            with open(fasta_file, encoding="utf-8") as f:
                return list(SeqIO.parse(f, "fasta"))
