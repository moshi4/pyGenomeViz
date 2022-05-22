from __future__ import annotations

from dataclasses import dataclass
from typing import List


class Track:
    def __init__(self, name: str, size: int, type: str = "feature"):
        """Constructor"""
        self.name = name
        self.size = size
        self.type = type
        self.features: List[Feature] = []

    def add_feature(
        self,
        start: int,
        end: int,
        strand: int,
        label: str = "",
        plotstyle: str = "bigarrow",
        facecolor: str = "orange",
        edgecolor: str = "black",
    ) -> None:
        """Add feature"""
        self.features.append(
            Feature(start, end, strand, label, plotstyle, facecolor, edgecolor)
        )


@dataclass
class Feature:
    """Feature DataClass"""

    start: int
    end: int
    strand: int  # -1 or 1
    label: str = ""
    plotstyle: str = "bigarrow"  # "bigarrow", "arrow", "box"
    facecolor: str = "orange"
    edgecolor: str = "black"

    def __post_init__(self):
        # Check start, end postion
        if self.start > self.end:
            err_msg = f"Feature 'end' must be larger than 'start' ({self})"
            raise ValueError(err_msg)
        # Check strand
        if self.strand not in (1, -1):
            err_msg = f"Strand must be '1' or '-1' ('{self.strand}')"
            raise ValueError(err_msg)
        # Check feature plot style
        self.plotstyle = self.plotstyle.lower()
        if self.plotstyle not in ("bigarrow", "arrow", "box"):
            err_msg = (
                f"'Style must be 'bigarrow' or 'arrow' or 'box' ('{self.plotstyle}')"
            )
            raise ValueError(err_msg)

    @property
    def length(self) -> int:
        """Feature length"""
        return self.end - self.start

    def __add__(self, offset: int) -> Feature:
        """Add offset"""
        return Feature(
            self.start + offset,
            self.end + offset,
            self.strand,
            self.label,
            self.plotstyle,
            self.facecolor,
            self.edgecolor,
        )
