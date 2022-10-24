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
            GFF file (`.gz`, `.bz2`, `.zip` format file is automatically uncompressed)
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
        self.target_seqid = target_seqid
        self._records = self._parse_gff(gff_file, min_range, max_range)

    def _parse_gff(
        self,
        gff_file: str | Path,
        min_range: int | None,
        max_range: int | None,
    ) -> list[GffRecord]:
        """Parse GFF

        Only parse target seqid record.
        If target_record is None, only parse first seqid record.

        Parameters
        ----------
        gff_file : str | Path
            GFF file
        min_range : int | None
            Min range
        max_range : int | None
            Max range

        Returns
        -------
        list[GffRecord]
            GFF record list
        """
        gff_file = Path(gff_file)
        if gff_file.suffix == ".gz":
            with gzip.open(gff_file, mode="rt") as f:
                gff_records, start, end = GffRecord._parse(f, self.target_seqid)
        elif gff_file.suffix == ".bz2":
            with bz2.open(gff_file, mode="rt") as f:
                gff_records, start, end = GffRecord._parse(f, self.target_seqid)
        elif gff_file.suffix == ".zip":
            with zipfile.ZipFile(gff_file) as zip:
                with zip.open(zip.namelist()[0]) as f:
                    io = TextIOWrapper(f)
                    gff_records, start, end = GffRecord._parse(io, self.target_seqid)
        else:
            with open(gff_file) as f:
                gff_records, start, end = GffRecord._parse(f, self.target_seqid)

        self.min_range = start if min_range is None else min_range
        self.max_range = end if max_range is None else max_range

        return gff_records

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
    def range_size(self) -> int:
        """Range size"""
        return self.max_range - self.min_range

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
            feature_props = dict(
                type=feature_type,
                strand=parent_record.strand,
                id=parent_record.attrs.get("ID", [""])[0],
                qualifiers=parent_record.attrs,
            )
            if len(exons) == 1:
                loc = FeatureLocation(exons[0].start, exons[0].end, exons[0].strand)
                exon_feature = SeqFeature(loc, **feature_props)
            elif len(exons) >= 2:
                exons = sorted(exons, key=lambda e: e.start)
                locs = [FeatureLocation(e.start, e.end, e.strand) for e in exons]
                exon_feature = SeqFeature(CompoundLocation(locs), **feature_props)
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
    def _parse(
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

        def is_gff_record_line(line: str) -> bool:
            if line.startswith("#") or len(line.split("\t")) < 9:
                return False
            else:
                return True

        gff_all_lines = handle.read().splitlines()
        gff_record_lines = [ln for ln in gff_all_lines if is_gff_record_line(ln)]
        gff_records = [GffRecord._parse_line(ln) for ln in gff_record_lines]

        seqid_list = list(dict.fromkeys([rec.seqid for rec in gff_records]))
        target_seqid = seqid_list[0] if target_seqid is None else target_seqid
        gff_records = [rec for rec in gff_records if rec.seqid == target_seqid]

        try:
            target = f"##sequence-region\t{target_seqid}\t"
            region_line = [ln for ln in gff_all_lines if ln.startswith(target)][0]
            region_elms = region_line.split("\t")
            start, end = int(region_elms[2]) - 1, int(region_elms[3])
        except Exception:
            start, end = 0, max([r.end for r in gff_records])

        return gff_records, start, end

    @staticmethod
    def _parse_line(gff_line: str) -> GffRecord:
        """Parse GFF record line

        Parameters
        ----------
        gff_line : str
            GFF record line

        Returns
        -------
        GffRecord
            Gff record (0-based index)
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
