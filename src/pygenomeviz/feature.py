from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Tuple


@dataclass
class Feature:
    """Feature DataClass"""

    start: int
    end: int
    strand: int  # -1 or 1
    label: str = ""
    labelsize: int = 15
    labelcolor: str = "black"
    plotstyle: str = "bigarrow"  # "bigarrow", "arrow", "bigbox", "box"
    facecolor: str = "orange"
    edgecolor: str = "black"
    linewidth: float = 0
    labelrotation: int = 0

    def __post_init__(self):
        # Change unknown strand value to 1
        if self.strand not in (1, -1):
            self.strand = 1
        # Check start, end postion
        if self.start > self.end:
            err_msg = f"Feature 'end' must be larger than 'start' ({self})"
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
        feature_size_ratio: float,
        arrow_shaft_ratio: float,
    ) -> Dict[str, Any]:
        """Feature object drawing parameters"""
        ylim = (ylim[0] * feature_size_ratio, ylim[1] * feature_size_ratio)
        # x, y
        x = self.end if self.strand == -1 else self.start
        if self.plotstyle in ("bigarrow", "bigbox"):
            y = 0
        else:
            if self.strand == -1:
                y = ylim[0] / 2
            else:
                y = ylim[1] / 2
        # dx, dy
        dx, dy = self.length * self.strand, 0
        # head width
        max_width = ylim[1] - ylim[0]  # = 2.0
        if self.plotstyle in ("bigarrow", "bigbox"):
            head_width = max_width
        else:
            head_width = max_width / 2
        # shaft_width
        if self.plotstyle in ("bigarrow", "arrow"):
            shaft_width = head_width * arrow_shaft_ratio
        else:
            shaft_width = head_width
        # head length
        if self.plotstyle in ("bigbox", "box"):
            head_length = 0
        else:
            head_length = max_track_size * 0.015
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
            "linewidth": self.linewidth,
        }

    def text_params(
        self, ylim: Tuple[float, float], feature_size_ratio: float
    ) -> Dict[str, Any]:
        """Feature text drawing parameters"""
        x = (self.start + self.end) / 2
        ylim = (ylim[0] * feature_size_ratio, ylim[1] * feature_size_ratio)
        # if self.plotstyle in ("bigarrow", "bigbox"):
        #     y = 0
        if self.strand == -1:
            # y = ylim[0] / 2
            y = ylim[0]
            ha, va = "left", "top"
            labelrotation = self.labelrotation * self.strand
        else:
            # y = ylim[1] / 2
            y = ylim[1]
            ha, va = "left", "bottom"
            labelrotation = self.labelrotation
        return {
            "x": x,
            "y": y,
            "s": self.label,
            "color": self.labelcolor,
            "fontsize": self.labelsize,
            "rotation": labelrotation,
            "ha": ha,
            "va": va,
            "zorder": 10,
            "rotation_mode": "anchor",  # 'anchor' or 'default'
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
            self.linewidth,
            self.labelrotation,
        )
