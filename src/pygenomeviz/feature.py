from __future__ import annotations

from dataclasses import dataclass


@dataclass
class Feature:
    """Feature DataClass"""

    start: int
    end: int
    strand: int  # -1 or 1
    label: str = ""
    labelsize: int = 15
    labelcolor: str = "black"
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
        return Feature(
            self.start + offset,
            self.end + offset,
            self.strand,
            self.label,
            self.labelsize,
            self.labelcolor,
            self.plotstyle,
            self.facecolor,
            self.edgecolor,
        )
