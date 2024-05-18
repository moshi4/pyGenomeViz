from __future__ import annotations

import csv
import io
from dataclasses import astuple, dataclass
from pathlib import Path

from pygenomeviz.typing import SeqType


@dataclass
class AlignCoord:
    """Alignment Coordinates DataClass (0-based coordinate)"""

    query_id: str
    query_name: str
    query_start: int
    query_end: int
    ref_id: str
    ref_name: str
    ref_start: int
    ref_end: int
    identity: float | None = None
    evalue: float | None = None

    @property
    def query_length(self) -> int:
        """Query length"""
        return abs(self.query_end - self.query_start)

    @property
    def query_strand(self) -> int:
        """Query strand"""
        return 1 if self.query_end >= self.query_start else -1

    @property
    def query_link(self) -> tuple[str, str, int, int]:
        """Query (name, start, end) link"""
        return (self.query_id, self.query_name, self.query_start, self.query_end)

    @property
    def query_block(self) -> tuple[int, int, int]:
        """Query (start, end, strand) block"""
        if self.query_start < self.query_end:
            return (self.query_start, self.query_end, self.query_strand)
        else:
            return (self.query_end, self.query_start, self.query_strand)

    @property
    def ref_length(self) -> int:
        """Reference length"""
        return abs(self.query_end - self.query_start)

    @property
    def ref_strand(self) -> int:
        """Reference strand"""
        return 1 if self.ref_end >= self.ref_start else -1

    @property
    def ref_link(self) -> tuple[str, str, int, int]:
        """Reference (name, start, end) link"""
        return (self.ref_id, self.ref_name, self.ref_start, self.ref_end)

    @property
    def ref_block(self) -> tuple[int, int, int]:
        """Reference (start, end, strand) block"""
        if self.ref_start < self.ref_end:
            return (self.ref_start, self.ref_end, self.ref_strand)
        else:
            return (self.ref_end, self.ref_start, self.ref_strand)

    @property
    def is_inverted(self) -> bool:
        """Check inverted or not"""
        return self.query_strand * self.ref_strand < 0

    @property
    def as_tsv_format(self) -> str:
        """TSV format text"""
        return "\t".join(
            (
                str(self.query_id),
                str(self.query_name),
                str(self.query_start),
                str(self.query_end),
                str(self.query_length),
                str(self.ref_id),
                str(self.ref_name),
                str(self.ref_start),
                str(self.ref_end),
                str(self.ref_length),
                str(self.identity) if self.identity is not None else "na",
                str(self.evalue) if self.evalue is not None else "na",
            )
        )

    @staticmethod
    def parse_blast_file(
        blast_file: str | Path,
        query_id: str,
        ref_id: str,
    ) -> list[AlignCoord]:
        """Parse blast result file (outfmt=6)

        Parameters
        ----------
        blast_file : str | Path
            Blast result file
        query_id : str
            Query ID
        ref_id : str
            Reference ID

        Returns
        -------
        align_coords : list[AlignCoord]
            Alignment coords
        """
        align_coords = []
        with open(blast_file) as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                qseqid, sseqid = row[0], row[1]
                # Blast  pident: 0 <= pident <= 100
                # MMseqs pident: 0 <= pident <= 1.0
                pident = float(row[2])
                if 0 <= pident <= 1.0:
                    # Convert [0-1] to [0-100] pident (e.g. 0.95215 -> 95.21)
                    pident = int(float(pident) * 10000) / 100
                qstart, qend, sstart, send = map(int, row[6:10])
                qstart, sstart = qstart - 1, sstart - 1  # 1-based to 0-based coordinate
                evalue = float(row[10])

                align_coords.append(
                    AlignCoord(
                        query_id,
                        qseqid,
                        qstart,
                        qend,
                        ref_id,
                        sseqid,
                        sstart,
                        send,
                        pident,
                        evalue,
                    ),
                )
        return align_coords

    @staticmethod
    def parse_mummer_file(
        mummer_file: str | Path,
        query_id: str,
        ref_id: str,
        seqtype: SeqType,
    ) -> list[AlignCoord]:
        """Parse MUMmer(nucmer|promer) result file

        Parameters
        ----------
        coords_tsv_file : str | Path
            MUMmer align coords file
        query_id : str
            Query ID
        ref_id : str
            Reference ID
        seqtype : SeqType
            `nucleotide` or `protein`

        Returns
        -------
        align_coords : list[AlignCoord]
            Align coord list
        """
        align_coords = []
        with open(mummer_file) as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                # Check read file contents & extract required row values
                if seqtype == "nucleotide":
                    if len(row) != 9:
                        err_msg = f"Invalid nucmer coords file '{mummer_file}'!!"
                        raise ValueError(err_msg)
                elif seqtype == "protein":
                    if len(row) != 13:
                        err_msg = f"Invalid promer coords file '{mummer_file}'!!"
                        raise ValueError(err_msg)
                    row = row[0:7] + row[11:13]
                else:
                    raise ValueError(f"Invalid seqtype '{seqtype}'!!")

                # Convert to correct value type (1-based to 0-based coordinate)
                ref_start, ref_end = int(row[0]) - 1, int(row[1])
                query_start, query_end = int(row[2]) - 1, int(row[3])
                identity = float(row[6])
                ref_name, query_name = str(row[7]), str(row[8])

                align_coords.append(
                    AlignCoord(
                        query_id,
                        query_name,
                        query_start,
                        query_end,
                        ref_id,
                        ref_name,
                        ref_start,
                        ref_end,
                        identity,
                    ),
                )

        return align_coords

    @staticmethod
    def parse_pmauve_file(
        bbone_file: str | Path,
        names: list[str],
        refid: int = 0,
    ) -> list[AlignCoord]:
        """Parse progressiveMauve bbone file

        Parameters
        ----------
        bbone_file : str | Path
            progressiveMauve bbone format file
        names : list[str]
            Sequence names
        refid : int, optional
            Reference genome index

        Returns
        -------
        align_coords : list[AlignCoord]
            Align coord list
        """
        with open(bbone_file) as f:
            reader = csv.reader(f, delimiter="\t")
            header_row = next(reader)
            genome_num = int(len(header_row) / 2)
            rows = []
            for row in reader:
                row = [int(col) for col in row]
                ref_idx = refid * 2
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
            rows = sorted(rows, key=lambda row: row[ref_idx])  # type: ignore

        align_coords = []
        for row in rows:
            for i in range(genome_num - 1):
                idx = i * 2
                # Query start-end
                qstart, qend = row[idx], row[idx + 1]
                if qstart < 0 and qend < 0:
                    qstart, qend = abs(qend), abs(qstart)
                # Reference start-end
                rstart, rend = row[idx + 2], row[idx + 3]
                if rstart < 0 and rend < 0:
                    rstart, rend = abs(rend), abs(rstart)
                qname, rname = names[i], names[i + 1]
                align_coord = AlignCoord(
                    qname,
                    qname,
                    qstart,
                    qend,
                    rname,
                    rname,
                    rstart,
                    rend,
                )
                align_coords.append(align_coord)
        return align_coords

    @staticmethod
    def write(
        align_coords: list[AlignCoord],
        outfile: str | Path | io.StringIO | io.BytesIO,
    ) -> None:
        """Write alignment coords as tsv format file

        Parameters
        ----------
        align_coords : list[AlignCoord]
            Alignment coords
        outfile : str | Path | StringIO | BytesIO
            Output file path or io stream
        """
        header_list = [
            "QUERY_ID",
            "QUERY_NAME",
            "QUERY_START",
            "QUERY_END",
            "QUERY_LENGTH",
            "REF_ID",
            "REF_NAME",
            "REF_START",
            "REF_END",
            "REF_LENGTH",
            "IDENTITY",
            "EVALUE",
        ]
        header = "\t".join(header_list)
        output = "\n".join([ac.as_tsv_format for ac in align_coords])
        contents = f"{header}\n{output}"
        if isinstance(outfile, io.StringIO):
            outfile.write(contents)
        elif isinstance(outfile, io.BytesIO):
            outfile.write(bytes(contents, encoding="utf-8"))
        else:
            with open(outfile, "w") as f:
                f.write(contents)

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
                    if idx in (0, 1, 5, 6):
                        # qid, qname, rid, rname
                        typed_row.append(str(val))
                    elif idx in (2, 3, 7, 8):
                        # qstart, qend, rstart, rend
                        typed_row.append(int(val))
                    elif idx in (10, 11):
                        # identity, evalue
                        typed_row.append(float(val) if val != "na" else None)
                align_coords.append(AlignCoord(*typed_row))
        return align_coords

    @staticmethod
    def filter(
        align_coords: list[AlignCoord],
        *,
        length_thr: int | None = None,
        identity_thr: float | None = None,
        evalue_thr: float | None = None,
    ) -> list[AlignCoord]:
        """Filter align coords by `length` & `identity` & `evalue`

        Parameters
        ----------
        align_coords : list[AlignCoord]
            Align coord list
        length_thr : int | None, optional
            Length filter threshold
        identity_thr : float | None, optional
            Identity filter threshold
        evalue_thr : float | None, optional
            E-value filter threshold

        Returns
        -------
        filtered_align_coords : list[AlignCoord]
            Filtered align coord list
        """
        filtered_align_coords: list[AlignCoord] = []
        for ac in align_coords:
            qlen, rlen = ac.query_length, ac.ref_length
            if length_thr and (qlen < length_thr or rlen < length_thr):
                continue
            if identity_thr and ac.identity and ac.identity < identity_thr:
                continue
            if evalue_thr and ac.evalue and ac.evalue > evalue_thr:
                continue
            filtered_align_coords.append(AlignCoord(*astuple(ac)))
        return filtered_align_coords

    def __contains__(self, target_ac: AlignCoord) -> bool:
        """Check whether target is completely overlapping with self"""
        ac1, ac2 = target_ac, self
        ac1_qmin = min(ac1.query_start, ac1.query_end)
        ac1_qmax = max(ac1.query_start, ac1.query_end)
        ac2_qmin = min(ac2.query_start, ac2.query_end)
        ac2_qmax = max(ac2.query_start, ac2.query_end)
        ac1_rmin = min(ac1.ref_start, ac1.ref_end)
        ac1_rmax = max(ac1.ref_start, ac1.ref_end)
        ac2_rmin = min(ac2.ref_start, ac2.ref_end)
        ac2_rmax = max(ac2.ref_start, ac2.ref_end)
        if (
            ac2_qmin <= ac1_qmin <= ac1_qmax <= ac2_qmax
            and ac2_rmin <= ac1_rmin <= ac1_rmax <= ac2_rmax
        ):
            return True
        else:
            return False
