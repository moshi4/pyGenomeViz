from functools import cached_property
from pathlib import Path
from typing import Any, List, Optional, Union

from Bio import SeqIO, SeqUtils
from Bio.SeqFeature import FeatureLocation, Seq, SeqFeature
from Bio.SeqRecord import SeqRecord


class Genbank:
    """Genbank Class"""

    def __init__(
        self,
        gbk_file: Union[str, Path],
        name: str = "",
    ):
        """
        Parameters
        ----------
        gbk_file : Union[str, Path]
            Genbank file
        name : str, optional
            name
        """
        self.gbk_file: Path = Path(gbk_file)
        self.name: str = name if name != "" else self.gbk_file.with_suffix("").name
        self._records: List[SeqRecord] = list(SeqIO.parse(gbk_file, "genbank"))

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

    def gc_skew(self, window_size: int = 5000, step_size: int = 2000) -> List[float]:
        """Calculate GC skew in sliding window

        Parameters
        ----------
        window_size : int, optional
            Window size
        step_size : int, optional
            Step size

        Returns
        -------
        gc_skews : List[float]
            GC skew value list
        """
        gc_skew_values = []
        seq = self.genome_seq
        for i in range(0, len(seq), step_size):
            start_pos = i - int(window_size / 2)
            start_pos = 0 if start_pos < 0 else start_pos
            end_pos = i + int(window_size / 2)
            end_pos = len(seq) if end_pos > len(seq) else end_pos

            subseq = seq[start_pos:end_pos]
            g = subseq.count("G") + subseq.count("g")
            c = subseq.count("C") + subseq.count("c")
            try:
                skew = (g - c) / float(g + c)
            except ZeroDivisionError:
                skew = 0.0
            gc_skew_values.append(skew)
        return gc_skew_values

    def gc_content(self, window_size: int = 5000, step_size: int = 2000) -> List[float]:
        """Calculate GC content in sliding window

        Parameters
        ----------
        window_size : int, optional
            Window size
        step_size : int, optional
            Step size

        Returns
        -------
        gc_contents : List[float]
            GC content value list
        """
        gc_content_values = []
        seq = self.genome_seq
        for i in range(0, len(seq), step_size):
            start_pos = i - int(window_size / 2)
            start_pos = 0 if start_pos < 0 else start_pos
            end_pos = i + int(window_size / 2)
            end_pos = len(seq) if end_pos > len(seq) else end_pos

            subseq = seq[start_pos:end_pos]
            gc_content_values.append(SeqUtils.GC(subseq))
        return gc_content_values

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
