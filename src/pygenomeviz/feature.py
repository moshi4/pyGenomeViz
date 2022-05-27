from __future__ import annotations

from copy import deepcopy
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
    labelrotation: int = 30
    labelvpos: str = "strand"  # "top", "center", "bottom", "strand"
    labelhpos: str = "center"  # "left", "center", "right"
    labelha: str = "left"  # "left", "center", "right"

    def __post_init__(self):
        # Change unknown strand value to 1
        if self.strand not in (1, -1):
            self.strand = 1
        # Check start, end postion
        if self.start > self.end:
            err_msg = f"Feature 'end' must be larger than 'start' ({self})"
            raise ValueError(err_msg)
        # Check labelvpos
        if self.labelvpos not in ("top", "center", "bottom", "strand"):
            raise ValueError(f"'labelvpos={self.labelvpos}' is invalid parameter.")
        # Check labelhpos
        if self.labelhpos not in ("left", "center", "right"):
            raise ValueError(f"'labelhpos={self.labelhpos}' is invalid parameter.")
        # Check labelha
        if self.labelha not in ("left", "center", "right"):
            raise ValueError(f"'labelha={self.labelha}' is invalid parameter.")
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
        max_width = ylim[1] - ylim[0]
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
        # x
        if self.labelhpos == "left":
            x = self.start
        elif self.labelhpos == "right":
            x = self.end
        else:  # "center"
            x = (self.start + self.end) / 2

        # labelrotation
        labelrotation = self.labelrotation
        if self.labelvpos == "strand":
            labelrotation = labelrotation * self.strand

        # labelvpos
        labelvpos = self.labelvpos
        if labelvpos == "strand":
            labelvpos = "bottom" if self.strand == -1 else "top"

        # labelva, y
        ylim = (ylim[0] * feature_size_ratio, ylim[1] * feature_size_ratio)
        if labelvpos == "top":
            labelva = "bottom"
            y = ylim[1]
        elif labelvpos == "bottom":
            labelva = "top"
            y = ylim[0]
        else:  # "center"
            labelva = "center"
            if self.plotstyle in ("bigarrow", "bigbox"):
                y = 0
            else:
                y = (abs(ylim[0]) / 2) * self.strand

        return {
            "x": x,
            "y": y,
            "s": self.label,
            "color": self.labelcolor,
            "fontsize": self.labelsize,
            "rotation": labelrotation,
            "ha": self.labelha,
            "va": labelva,
            "zorder": 10,
            "rotation_mode": "anchor",
        }

    def __add__(self, offset: int):
        feature = deepcopy(self)
        feature.start += offset
        feature.end += offset
        return feature
