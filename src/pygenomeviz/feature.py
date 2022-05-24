from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Tuple


@dataclass
class Feature:
    """Feature DataClass"""

    start: int
    end: int
    strand: int  # -1 or 0 or 1
    label: str = ""
    labelsize: int = 15
    labelcolor: str = "black"
    plotstyle: str = "bigarrow"  # "bigarrow", "arrow", "bigbox", "box"
    facecolor: str = "orange"
    edgecolor: str = "black"
    labelrotation: int = 0

    def __post_init__(self):
        # Check start, end postion
        if self.start > self.end:
            err_msg = f"Feature 'end' must be larger than 'start' ({self})"
            raise ValueError(err_msg)
        # Check strand
        if self.strand not in (1, 0, -1):
            err_msg = f"Strand must be '1 | 0 | -1' (strand={self.strand})"
            raise ValueError(err_msg)
        # Check feature plot style
        self.plotstyle = self.plotstyle.lower()
        if self.plotstyle not in ("bigarrow", "arrow", "bigbox", "box"):
            err_msg = f"'Style must be 'bigarrow|arrow|bigbox|box' ('{self.plotstyle}')"
            raise ValueError(err_msg)

    @property
    def length(self) -> int:
        """Feature length"""
        return self.end - self.start

    def obj_params(
        self,
        max_track_size: int,
        ylim: Tuple[float, float],
        feature_track_pad: float,
        arrow_shaft_ratio: float,
    ) -> Dict[str, Any]:
        """Feature object drawing parameters"""
        # x, y
        x = self.end if self.strand == -1 else self.start
        if self.plotstyle in ("bigarrow", "bigbox"):
            y = 0
        else:
            if self.strand == -1:
                y = -0.5 + (feature_track_pad / 2)
            else:
                y = 0.5 - (feature_track_pad / 2)
        # dx, dy
        dx = -self.length if self.strand == -1 else self.length
        dy = 0
        # head width
        max_width = ylim[1] - ylim[0]  # = 2.0
        if self.plotstyle in ("bigarrow", "bigbox"):
            head_width = max_width - (feature_track_pad * 2)
        else:
            head_width = (max_width / 2) - feature_track_pad
        # shaft_width
        if self.plotstyle in ("bigarrow", "arrow"):
            shaft_width = head_width * arrow_shaft_ratio
        else:
            shaft_width = head_width
        # head length
        if self.plotstyle in ("bigbox", "box"):
            head_length = 0
        else:
            head_length = max_track_size * 0.02
            if abs(self.length) < head_length:
                head_length = abs(self.length)
        # zorder
        zorder = 5 if self.plotstyle in ("bigarrow", "bigbox") else -5

        return {
            "x": x,
            "y": y,
            "dx": dx,
            "dy": dy,
            "width": shaft_width,
            "head_width": head_width,
            "head_length": head_length,
            "fc": self.facecolor,
            "ec": self.edgecolor,
            "length_includes_head": True,
            "zorder": zorder,
        }

    @property
    def text_params(self) -> Dict[str, Any]:
        """Feature text drawing parameters"""
        return {
            "s": self.label,
            "color": self.labelcolor,
            "fontsize": self.labelsize,
            "rotation": self.labelrotation,
            "ha": "center",
            "va": "center",
            "zorder": 10,
        }

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
            self.labelrotation,
        )
