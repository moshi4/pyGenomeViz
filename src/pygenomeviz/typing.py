from __future__ import annotations

from typing import Literal

TrackAlignType = Literal["left", "center", "right"]
PlotStyle = Literal["bigarrow", "arrow", "bigbox", "box", "bigrbox", "rbox"]
VPos = Literal["top", "center", "bottom"]
HPos = Literal["left", "center", "right"]
SeqType = Literal["nucleotide", "protein"]
Unit = Literal["Gb", "Mb", "Kb", "bp"]

AlnCliName = Literal["pgv-blast", "pgv-mummer", "pgv-mmseqs", "pgv-pmauve"]
GenbankDatasetName = Literal[
    "acinetobacter_phage",
    "yersinia_phage",
    "enterobacteria_phage",
    "mycoplasma_mycoides",
    "escherichia_coli",
    "saccharomyces",
]
GffExampleFileName = Literal[
    "enterobacteria_phage.gff",
    "mycoplasma_mycoides.gff",
    "escherichia_coli.gff.gz",
    "saccharomyces_cerevisiae.gff.gz",
]

# GUI
AlnMethod = Literal[
    "MUMmer (nucleotide)",
    "MUMmer (protein)",
    "MMseqs RBH",
    "BLAST (nucleotide)",
    "BLAST (protein)",
    None,
]
