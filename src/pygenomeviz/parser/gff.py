from __future__ import annotations

import bz2
import gzip
import zipfile
from collections import defaultdict
from dataclasses import dataclass
from io import TextIOWrapper
from pathlib import Path
from typing import Any, TextIO

from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature


class Gff:
    """GFF Class"""

    def __init__(
        self,
        gff_file: str | Path,
        name: str | None = None,
        target_seqid: str | None = None,
        min_range: int | None = None,
        max_range: int | None = None,
    ):
        """
        Parameters
        ----------
        gff_file : str | Path
            GFF file (`.gz`, `.bz2`, `.zip` compressed file can be readable)
        name : str | None, optional
            name (If None, `file name` is set)
        target_seqid : str | None, optional
            Target seqid to be extracted. If None, only first seqid record is extracted.
        min_range : int | None, optional
            Min range to be extracted
        max_range : int | None, optional
            Max range to be extracted
        """
        self._gff_file = Path(gff_file)
        self._name = name
        self._target_seqid = target_seqid
        self._records, start, end = self._parse_gff(gff_file, target_seqid)
        self.min_range = start if min_range is None else min_range
        self.max_range = end if max_range is None else max_range
        if not 0 <= self.min_range <= self.max_range:
            err_msg = "Range must be '0 <= min_range <= max_range' "
            err_msg += f"(min_range={self.min_range}, max_range={self.max_range})"
            raise ValueError(err_msg)

    def _parse_gff(
        self,
        gff_file: str | Path,
        target_seqid: str | None,
    ) -> tuple[list[GffRecord], int, int]:
        """Parse GFF file

        Only parse target seqid record.
        If target_record is None, only parse first seqid record.

        Parameters
        ----------
        gff_file : str | Path
            GFF file
        target_seqid : str | None
            Target seqid to be extracted

        Returns
        -------
        tuple[list[GffRecord], int, int]
            GFF record list, start, end
        """
        gff_file = Path(gff_file)
        if gff_file.suffix == ".gz":
            with gzip.open(gff_file, mode="rt") as f:
                gff_records, start, end = self._parse_gff_textio(f, target_seqid)
        elif gff_file.suffix == ".bz2":
            with bz2.open(gff_file, mode="rt") as f:
                gff_records, start, end = self._parse_gff_textio(f, target_seqid)
        elif gff_file.suffix == ".zip":
            with zipfile.ZipFile(gff_file) as zip:
                with zip.open(zip.namelist()[0]) as f:
                    io = TextIOWrapper(f)
                    gff_records, start, end = self._parse_gff_textio(io, target_seqid)
        else:
            with open(gff_file) as f:
                gff_records, start, end = self._parse_gff_textio(f, target_seqid)

        return gff_records, start, end

    def _parse_gff_textio(
        self,
        handle: TextIO,
        target_seqid: str | None = None,
    ) -> tuple[list[GffRecord], int, int]:
        """Parse GFF file TextIO

        Parameters
        ----------
        handle : TextIO
            GFF TextIO handle
        target_seqid : str | None, optional
            GFF target seqid

        Returns
        -------
        tuple[list[GffRecord], int, int]
            GFF record list, start, end
        """
        gff_all_lines = handle.read().splitlines()
        gff_record_lines = filter(GffRecord.is_gff_record_line, gff_all_lines)
        gff_records = list(map(GffRecord.parse_line, gff_record_lines))
        if len(gff_records) == 0:
            err_msg = f"Failed to parse '{self._gff_file}' as GFF file "
            raise ValueError(err_msg)

        seqid_list = list(dict.fromkeys([rec.seqid for rec in gff_records]))
        self._seqid_list = seqid_list
        if target_seqid is None:
            target_seqid = seqid_list[0]
        if target_seqid not in seqid_list:
            err_msg = f"Not found target_seqid='{target_seqid}' in '{self._gff_file}'"
            raise ValueError(err_msg)
        gff_records = [rec for rec in gff_records if rec.seqid == target_seqid]

        try:
            target = f"##sequence-region\t{target_seqid}\t"
            region_line = [ln for ln in gff_all_lines if ln.startswith(target)][0]
            region_elms = region_line.split("\t")
            start, end = int(region_elms[2]) - 1, int(region_elms[3])
        except Exception:
            start, end = 0, max([r.end for r in gff_records])

        return gff_records, start, end

    @property
    def name(self) -> str:
        """Name"""
        if self._name is not None:
            return self._name
        if self._gff_file.suffix in (".gz", ".bz2", ".zip"):
            return self._gff_file.with_suffix("").with_suffix("").name
        else:
            return self._gff_file.with_suffix("").name

    @property
    def records(self) -> list[GffRecord]:
        """GFF records"""
        return self._records

    @property
    def range_size(self) -> int:
        """Range size (`max_range - min_range`)"""
        return self.max_range - self.min_range

    @property
    def seqid_list(self) -> list[str]:
        """seqid list"""
        return self._seqid_list

    def extract_features(self, feature_type: str = "CDS") -> list[SeqFeature]:
        """Extract features

        Parameters
        ----------
        feature_type : str, optional
            Feature type (`CDS`, `gene`, `mRNA`, etc...)

        Returns
        -------
        list[SeqFeature]
            Feature list
        """
        gff_records = [rec for rec in self._records if rec.type == feature_type]
        features: list[SeqFeature] = []
        for rec in gff_records:
            if not rec.is_within_range(self.min_range, self.max_range):
                continue
            feature = SeqFeature(
                FeatureLocation(rec.start, rec.end, rec.strand),
                type=rec.type,
                strand=rec.strand,
                id=rec.attrs.get("ID", [""])[0],
                qualifiers=rec.attrs,
            )
            features.append(feature)

        return features

    def extract_exon_features(self, feature_type: str = "mRNA") -> list[SeqFeature]:
        """Extract exon structure features

        Extract exons based on `parent feature` and `exon` ID-Parent relation

        Parameters
        ----------
        feature_type : str, optional
            Feature type (e.g. `mRNA`, `ncRNA` , etc...)

        Returns
        -------
        list[SeqFeature]
            Feature list
        """
        # Extract exon features by mRNA-exon relation
        parent_id = None
        parent_id2record: dict[str, GffRecord] = {}
        parent_id2exons: dict[str, list[GffRecord]] = defaultdict(list)
        for rec in self._records:
            if not rec.is_within_range(self.min_range, self.max_range):
                continue
            if rec.type == feature_type:
                parent_id = rec.attrs["ID"][0]
                parent_id2record[parent_id] = rec
            if rec.type == "exon":
                if parent_id is not None and parent_id == rec.attrs["Parent"][0]:
                    parent_id2exons[parent_id].append(rec)

        # Set exon features
        exon_features: list[SeqFeature] = []
        for parent_id in parent_id2record.keys():
            parent_record = parent_id2record[parent_id]
            exons = parent_id2exons[parent_id]
            feature_kws = dict(
                type=feature_type,
                strand=parent_record.strand,
                id=parent_record.attrs.get("ID", [""])[0],
                qualifiers=parent_record.attrs,
            )
            if len(exons) == 1:
                loc = FeatureLocation(exons[0].start, exons[0].end, exons[0].strand)
                exon_feature = SeqFeature(loc, **feature_kws)
            elif len(exons) >= 2:
                exons = sorted(exons, key=lambda e: e.start)
                locs = [FeatureLocation(e.start, e.end, e.strand) for e in exons]
                exon_feature = SeqFeature(CompoundLocation(locs), **feature_kws)
            else:
                # If no exon exists, skip feature extraction
                continue

            exon_features.append(exon_feature)

        return exon_features

    def __str__(self):
        return f"{self.name} ({self.min_range:,} - {self.max_range:,} bp)"


@dataclass
class GffRecord:
    seqid: str
    source: str
    type: str
    start: int
    end: int
    score: float | None
    strand: int
    frame: int | None
    attrs: dict[str, list[str]]

    def is_within_range(self, min_range: int, max_range: int) -> bool:
        """Check within target range or not

        Parameters
        ----------
        min_range : int
            Min range
        max_range : int
            Max range

        Returns
        -------
        bool
            Check result
        """
        if min_range <= self.start <= self.end <= max_range:
            return True
        else:
            return False

    @staticmethod
    def is_gff_record_line(line: str) -> bool:
        """Check GFF record line or not

        Parameters
        ----------
        line : str
            GFF line

        Returns
        -------
        bool
            Check result
        """
        if line.startswith("#") or len(line.split("\t")) < 9:
            return False
        else:
            return True

    @staticmethod
    def parse_line(gff_line: str) -> GffRecord:
        """Parse GFF record line

        Parameters
        ----------
        gff_line : str
            GFF record line

        Returns
        -------
        GffRecord
            GFF record (0-based index)
        """
        gff_elms: list[Any] = gff_line.split("\t")[0:9]
        # start, end
        gff_elms[3], gff_elms[4] = int(gff_elms[3]) - 1, int(gff_elms[4])
        # score
        gff_elms[5] = None if gff_elms[5] == "." else float(gff_elms[5])
        # strand
        if gff_elms[6] == "+":
            gff_elms[6] = 1
        elif gff_elms[6] == "-":
            gff_elms[6] = -1
        else:
            gff_elms[6] = 0
        # frame
        gff_elms[7] = None if gff_elms[7] == "." else int(gff_elms[7])
        # attrs
        attr_dict: dict[str, list[str]] = {}
        attrs = str(gff_elms[8]).split(";")
        for attr in attrs:
            if attr != "":
                key, value = attr.split("=")
                attr_dict[key] = value.split(",")
        gff_elms[8] = attr_dict

        return GffRecord(*gff_elms)
