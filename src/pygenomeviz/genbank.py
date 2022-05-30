from functools import cached_property
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
    ):
        """
        Parameters
        ----------
        gbk_source : Union[str, Path, TextIOWrapper]
            Genbank file or source

        name : Optional[str]
            name (If None, `file name` or `record name` is set)
        """
        self._gbk_source = gbk_source
        self._name = name
        self._records: List[SeqRecord] = list(SeqIO.parse(gbk_source, "genbank"))

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
    def genome_length(self) -> int:
        """Genome sequence length"""
        return len(self.genome_seq)

    @property
    def genome_seq(self) -> str:
        """Genome sequence (join all contig sequences)"""
        return "".join(self.contig_seqs)

    @property
    def contig_seqs(self) -> List[str]:
        """Contig sequences"""
        return [str(r.seq) for r in self._records]

    @cached_property
    def average_gc(self) -> float:
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
        print(window_size, step_size)
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

    def extract_all_features(
        self,
        feature_type: str = "CDS",
        target_strand: Optional[int] = None,
    ) -> List[SeqFeature]:
        """Extract all features

        Parameters
        ----------
        feature_type : str, optional
            Extract feature type
        target_strand : Optional[int], optional
            Extract target strand

        Returns
        -------
        all_features : List[SeqFeature]
            Extracted all features
        """
        extract_features = []
        base_len = 0
        for record in self._records:
            features = [f for f in record.features if f.type == feature_type]
            for f in features:
                if feature_type == "CDS":
                    # Exclude pseudogene (no translated gene)
                    translation = f.qualifiers.get("translation", [None])[0]
                    if translation is None:
                        continue
                start = self._to_int(f.location.parts[0].start) + base_len
                end = self._to_int(f.location.parts[-1].end) + base_len
                # Exclude feature that straddle start position
                if start > end:
                    continue
                # Extract only target strand feature
                if target_strand is not None and f.strand != target_strand:
                    continue

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
        features = self.extract_all_features("CDS", None)
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
