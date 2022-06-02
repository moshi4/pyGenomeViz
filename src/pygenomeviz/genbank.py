from functools import lru_cache
from io import TextIOWrapper
from pathlib import Path
from typing import Any, List, Optional, Tuple, Union

import numpy as np
from Bio import SeqIO, SeqUtils
from Bio.SeqFeature import FeatureLocation, Seq, SeqFeature
from Bio.SeqRecord import SeqRecord


class Genbank:
    """Genbank Class"""

    def __init__(
        self,
        gbk_source: Union[str, Path, TextIOWrapper],
        name: Optional[str] = None,
        reverse: bool = False,
        min_range: Optional[int] = None,
        max_range: Optional[int] = None,
    ):
        """
        Parameters
        ----------
        gbk_source : Union[str, Path, TextIOWrapper]
            Genbank file or source
        name : Optional[str]
            name (If None, `file name` or `record name` is set)
        reverse : bool, optional
            Reverse genome or not
        min_range : Optional[int], optional
            Min range to be extracted (Default: `1`)
        max_range : Optional[int], optional
            Max range to be extracted (Default: `genome length`)
        """
        self._gbk_source = gbk_source
        self._name = name
        self._records: List[SeqRecord] = list(SeqIO.parse(gbk_source, "genbank"))
        self.reverse = reverse
        self.min_range = 1 if min_range is None else min_range
        self.max_range = self.full_genome_length if max_range is None else max_range

        if not 1 <= self.min_range <= self.max_range <= self.full_genome_length:
            err_msg = f"min_range={min_range}, max_range={max_range} is invalid. \n"
            err_msg += "Range must be "
            err_msg += f"'1 <= min_range <= max_range <= {self.full_genome_length}'"
            raise ValueError(err_msg)

    @property
    def name(self) -> str:
        """Name"""
        if self._name is not None:
            return self._name
        if isinstance(self._gbk_source, (str, Path)):
            return Path(self._gbk_source).with_suffix("").name
        elif isinstance(self._gbk_source, TextIOWrapper):
            return self._records[0].name
        else:
            raise NotImplementedError()

    @property
    def records(self) -> List[SeqRecord]:
        """Genbank records"""
        if self.reverse:
            return list(reversed([r.reverse_complement() for r in self._records]))
        else:
            return self._records

    @property
    def full_genome_length(self) -> int:
        """Full genome sequence length"""
        return len(self.full_genome_seq)

    @property
    def genome_length(self) -> int:
        """Genome sequence length"""
        return len(self.genome_seq)

    @property
    def full_genome_seq(self) -> str:
        """Genome sequence (join all contig sequences)"""
        return "".join(str(r.seq) for r in self.records)

    @property
    def genome_seq(self) -> str:
        """Genome sequence (join all contig sequences)"""
        seq = "".join(str(r.seq) for r in self.records)
        return seq[self.min_range - 1 : self.max_range]

    @lru_cache(maxsize=None)
    def calc_average_gc(self) -> float:
        """Average GC content"""
        return SeqUtils.GC(self.genome_seq)

    def calc_gc_skew(
        self,
        window_size: Optional[int] = None,
        step_size: Optional[int] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Calculate GC skew in sliding window

        Parameters
        ----------
        window_size : int, optional
            Window size (Default: `genome_size / 400`)
        step_size : int, optional
            Step size (Default: `genome_size / 1000`)

        Returns
        -------
        (pos_list, gc_skew_list) : Tuple[np.ndarray, np.ndarray]
            Position list & GC skew list
        """
        pos_list, gc_skew_list = [], []
        seq = self.genome_seq
        if window_size is None:
            window_size = int(len(seq) / 400)
        if step_size is None:
            step_size = int(len(seq) / 1000)
        for i in range(0, len(seq), step_size):
            start_pos = i - int(window_size / 2)
            start_pos = 0 if start_pos < 0 else start_pos
            end_pos = i + int(window_size / 2)
            end_pos = len(seq) if end_pos > len(seq) else end_pos
            middle_pos = int((end_pos + start_pos) / 2)
            pos_list.append(middle_pos)

            subseq = seq[start_pos:end_pos]
            g = subseq.count("G") + subseq.count("g")
            c = subseq.count("C") + subseq.count("c")
            try:
                skew = (g - c) / float(g + c)
            except ZeroDivisionError:
                skew = 0.0
            gc_skew_list.append(skew)

        return (np.array(pos_list), np.array(gc_skew_list))

    def calc_gc_content(
        self,
        window_size: Optional[int] = None,
        step_size: Optional[int] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Calculate GC content in sliding window

        Parameters
        ----------
        window_size : int, optional
            Window size (Default: `genome_size / 400`)
        step_size : int, optional
            Step size (Default: `genome_size / 1000`)

        Returns
        -------
        (pos_list, gc_content_list) : Tuple[np.ndarray, np.ndarray]
            Position list & GC content list
        """
        pos_list, gc_content_list = [], []
        seq = self.genome_seq
        if window_size is None:
            window_size = int(len(seq) / 400)
        if step_size is None:
            step_size = int(len(seq) / 1000)
        for i in range(0, len(seq), step_size):
            start_pos = i - int(window_size / 2)
            start_pos = 0 if start_pos < 0 else start_pos
            end_pos = i + int(window_size / 2)
            end_pos = len(seq) if end_pos > len(seq) else end_pos
            middle_pos = int((end_pos + start_pos) / 2)
            pos_list.append(middle_pos)

            subseq = seq[start_pos:end_pos]
            gc_content_list.append(SeqUtils.GC(subseq))

        return (np.array(pos_list), np.array(gc_content_list))

    def extract_features(
        self,
        feature_type: str = "CDS",
        target_strand: Optional[int] = None,
        fix_position: bool = True,
        partial_range: bool = False,
    ) -> List[SeqFeature]:
        """Extract features

        Parameters
        ----------
        feature_type : str, optional
            Extract feature type
        target_strand : Optional[int], optional
            Extract target strand
        fix_position : bool, optional
            Fix feature start & end position by specified min_range parameter
            (fixed_start = start - min_range - 1, fixed_end = end - min_range - 1)
        partial_range : bool, optional
            If True, features that are partially included in range are also extracted

        Returns
        -------
        features : List[SeqFeature]
            Extracted features
        """
        extract_features = []
        min_range, max_range = self.min_range - 1, self.max_range
        base_len = 0
        for record in self.records:
            features = [f for f in record.features if f.type == feature_type]
            for f in features:
                if feature_type == "CDS":
                    # Exclude pseudogene (no translated gene)
                    translation = f.qualifiers.get("translation", [None])[0]
                    if translation is None:
                        continue
                start = self._to_int(f.location.parts[0].start) + base_len
                end = self._to_int(f.location.parts[-1].end) + base_len
                # Restrict features in range
                if partial_range:
                    if (
                        not min_range <= start <= max_range
                        and not min_range <= end <= max_range
                    ):
                        continue
                else:
                    if not min_range <= start <= end <= max_range:
                        continue
                # Extract only target strand feature
                if target_strand is not None and f.strand != target_strand:
                    continue
                # Fix start & end position by min_range
                if fix_position:
                    start -= min_range
                    end -= min_range

                extract_features.append(
                    SeqFeature(
                        location=FeatureLocation(start, end, f.strand),
                        type=f.type,
                        qualifiers=f.qualifiers,
                    ),
                )
            base_len += len(record.seq)

        return extract_features

    def write_cds_fasta(
        self,
        fasta_outfile: Union[str, Path],
    ):
        """Write CDS protein features fasta file

        Parameters
        ----------
        fasta_outfile : Union[str, Path]
            CDS fasta file
        """
        features = self.extract_features("CDS", None, fix_position=False)
        cds_seq_records: List[SeqRecord] = []
        for idx, feature in enumerate(features, 1):
            qualifiers = feature.qualifiers
            protein_id = qualifiers.get("protein_id", [None])[0]
            product = qualifiers.get("product", [""])[0]
            translation = qualifiers.get("translation", [None])[0]

            start = self._to_int(feature.location.start)
            end = self._to_int(feature.location.end)
            strand = "-" if feature.strand == -1 else "+"

            location_id = f"|{start}_{end}_{strand}|"
            if protein_id is None:
                seq_id = f"GENE{idx:06d}{location_id}"
            else:
                seq_id = f"GENE{idx:06d}_{protein_id}{location_id}"

            cds_seq_record = SeqRecord(
                seq=Seq(translation), id=seq_id, description=product
            )
            cds_seq_records.append(cds_seq_record)

        SeqIO.write(cds_seq_records, fasta_outfile, "fasta-2line")

    def write_genome_fasta(
        self,
        outfile: Union[str, Path],
    ) -> None:
        """Write genome fasta file

        Parameters
        ----------
        outfile : Union[str, Path]
            Output genome fasta file
        """
        write_seq = self.genome_seq
        with open(outfile, "w") as f:
            f.write(f">{self.name}\n{write_seq}\n")

    def _to_int(self, value: Any) -> int:
        """Convert to int (Required for AbstractPostion|ExactPostion)"""
        return int(str(value).replace("<", "").replace(">", ""))
