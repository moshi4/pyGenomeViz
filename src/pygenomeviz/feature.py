from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

from matplotlib.patches import FancyArrow, Patch, PathPatch, Rectangle
from matplotlib.path import Path


@dataclass
class Feature:
    """Feature DataClass"""

    start: int
    end: int
    strand: int  # -1 or 1
    label: str = ""
    labelsize: int = 15
    labelcolor: str = "black"
    plotstyle: str = "bigarrow"  # "(big)arrow", "(big)box", "(big)rbox"
    facecolor: str = "orange"
    edgecolor: str = "black"
    linewidth: float = 0
    labelrotation: int = 45
    labelvpos: str = "strand"  # "top", "center", "bottom", "strand"
    labelhpos: str = "center"  # "left", "center", "right"
    labelha: str = "left"  # "left", "center", "right"
    arrow_shaft_ratio: float = 0.5
    size_ratio: float = 0.9
    exon_regions: Optional[List[Tuple[int, int]]] = None

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
        self.plotstyle_list = ("bigarrow", "arrow", "bigbox", "box", "bigrbox", "rbox")
        self.plotstyle = self.plotstyle.lower()
        if self.plotstyle not in self.plotstyle_list:
            err_msg = f"'plotstyle must be '{'|'.join(self.plotstyle_list)}'.\n"
            err_msg += f"plotstyle='{self.plotstyle}' is invalid."
            raise ValueError(err_msg)
        # Check arrow shaft ratio
        if not 0 <= self.arrow_shaft_ratio <= 1:
            err_msg = "'arrow_shaft_ratio' must be '0 <= value <= 1' "
            err_msg += f"(value={self.arrow_shaft_ratio})"
            raise ValueError(err_msg)
        # Check size ratio
        if not 0 <= self.size_ratio <= 1:
            err_msg = "'size_ratio' must be '0 <= value <= 1' "
            err_msg += f"(value={self.size_ratio})"
            raise ValueError(err_msg)
        # Check exon regions
        if self.exon_regions is not None:
            if len(self.exon_regions) < 2:
                raise ValueError()
            start_exon_loc = self.exon_regions[0][0]
            end_exon_loc = self.exon_regions[-1][-1]
            if self.start != start_exon_loc or self.end != end_exon_loc:
                raise ValueError()

    @property
    def length(self) -> int:
        """Feature length"""
        return self.end - self.start + 1

    @property
    def is_bigstyle(self) -> bool:
        """Check plotstyle is 'big~~~' or not"""
        return self.plotstyle.startswith("big")

    @property
    def patch_kwargs(self) -> Dict[str, Any]:
        """Patch keyword arguments dict"""
        return {
            "fc": self.facecolor,
            "ec": self.edgecolor,
            "lw": self.linewidth,
            "zorder": 5 if self.is_bigstyle else -5,
        }

    def box_patch(
        self, start: int, end: int, ylim: Tuple[float, float], ratio=1.0
    ) -> Rectangle:
        """Box patch"""
        # x, y
        x = start
        if self.is_bigstyle or self.strand == -1:
            y = ylim[0]
        else:
            y = 0
        y *= ratio
        # width, height
        width = end - start + 1
        if self.is_bigstyle:
            height = ylim[1] - ylim[0]
        else:
            height = ylim[1]
        height *= ratio

        return Rectangle((x, y), width, height, **self.patch_kwargs)

    def arrow_patch(
        self,
        start: int,
        end: int,
        ylim: Tuple[float, float],
        max_track_size: int,
        exon_box: bool = False,
    ) -> FancyArrow:
        """Arrow patch"""
        # x, y
        x = end if self.strand == -1 else start
        if self.is_bigstyle:
            y = 0
        else:
            if self.strand == -1:
                y = ylim[0] / 2
            else:
                y = ylim[1] / 2
        # dx, dy
        length = end - start + 1
        dx, dy = length * self.strand, 0
        # head width
        max_width = ylim[1] - ylim[0]
        head_width = max_width
        if self.is_bigstyle:
            head_width = max_width
        else:
            head_width = max_width / 2
        # shaft_width
        shaft_width = head_width * self.arrow_shaft_ratio
        # head length
        head_length = max_track_size * 0.015
        if abs(self.length) < head_length:
            head_length = abs(self.length)

        if exon_box:
            head_length = 0

        return FancyArrow(
            x,
            y,
            dx,
            dy,
            width=shaft_width,
            length_includes_head=True,
            head_width=head_width,
            head_length=head_length,
            **self.patch_kwargs,
        )

    def rbox_patch(
        self, start: int, end: int, ylim: Tuple[float, float], max_track_size: int
    ) -> PathPatch:
        """Rounded box patch"""
        length = end - start + 1
        r_size = max_track_size * 0.005
        if length <= 4 * r_size:
            r_size = length * 0.25

        xmin, xmax = start, end
        if self.is_bigstyle:
            ymin, ymax, ycenter = ylim[0], ylim[1], 0
        else:
            if self.strand == -1:
                ymin, ymax, ycenter = ylim[0], 0, ylim[0] / 2
            else:
                ymin, ymax, ycenter = 0, ylim[1], ylim[1] / 2

        path_data = [
            (Path.MOVETO, (xmin + r_size, ymax)),
            (Path.LINETO, (xmax - r_size, ymax)),
            (Path.CURVE3, (xmax + r_size, ycenter)),
            (Path.CURVE3, (xmax - r_size, ymin)),
            (Path.LINETO, (xmin + r_size, ymin)),
            (Path.CURVE3, (xmin - r_size, ycenter)),
            (Path.CURVE3, (xmin + r_size, ymax)),
        ]
        codes, verts = zip(*path_data)
        return PathPatch(Path(verts, codes), **self.patch_kwargs)

    def intron_patch(
        self, start: int, end: int, ylim: Tuple[float, float]
    ) -> PathPatch:
        """Intron patch"""
        xmin, xmax, xcenter = start, end, (start + end) / 2
        if self.is_bigstyle:
            ymin, ymax, ycenter = ylim[0], ylim[1], 0
        else:
            if self.strand == -1:
                ymin, ymax, ycenter = ylim[0], 0, ylim[0] / 2
            else:
                ymin, ymax, ycenter = 0, ylim[1], ylim[1] / 2
        ytop = ymin if self.strand == -1 else ymax

        path_data = [
            (Path.MOVETO, (xmin, ycenter)),
            (Path.LINETO, (xcenter, ytop)),
            (Path.LINETO, (xmax, ycenter)),
        ]
        codes, verts = zip(*path_data)
        return PathPatch(Path(verts, codes), lw=1, fill=False, zorder=5)

    @property
    def intron_regions(self) -> List[Tuple[int, int]]:
        """Intron regions"""
        intron_regions = []
        if self.exon_regions is None:
            return intron_regions
        for i in range(0, len(self.exon_regions) - 1):
            intron_start = self.exon_regions[i][1] + 1
            intron_end = self.exon_regions[i + 1][0] - 1
            intron_regions.append((intron_start, intron_end))
        return intron_regions

    def patches(
        self,
        max_track_size: int,
        ylim: Tuple[float, float],
    ) -> List[Patch]:
        """Feature patch"""
        ylim = (ylim[0] * self.size_ratio, ylim[1] * self.size_ratio)

        if self.exon_regions is not None:
            patches = []
            exon_regions = self.exon_regions[:: self.strand]
            for idx, (exon_start, exon_end) in enumerate(exon_regions, 1):
                if self.plotstyle in ("bigbox", "box"):
                    p = self.box_patch(exon_start, exon_end, ylim)
                elif self.plotstyle in ("bigarrow", "arrow"):
                    exon_box = False if idx == len(exon_regions) else True
                    p = self.arrow_patch(
                        exon_start, exon_end, ylim, max_track_size, exon_box
                    )
                elif self.plotstyle in ("bigrbox", "rbox"):
                    p = self.rbox_patch(exon_start, exon_end, ylim, max_track_size)
                else:
                    raise NotImplementedError()
                patches.append(p)

            for (intron_start, intron_end) in self.intron_regions:
                patches.append(self.intron_patch(intron_start, intron_end, ylim))

            return patches

        if self.plotstyle in ("bigbox", "box"):
            return [self.box_patch(self.start, self.end, ylim)]

        elif self.plotstyle in ("bigarrow", "arrow"):
            return [self.arrow_patch(self.start, self.end, ylim, max_track_size)]

        elif self.plotstyle in ("bigrbox", "rbox"):
            return [self.rbox_patch(self.start, self.end, ylim, max_track_size)]

        else:
            raise NotImplementedError()

    def text_params(self, ylim: Tuple[float, float]) -> Dict[str, Any]:
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
        ylim = (ylim[0] * self.size_ratio, ylim[1] * self.size_ratio)
        if labelvpos == "top":
            labelva = "bottom"
            y = ylim[1]
        elif labelvpos == "bottom":
            labelva = "top"
            y = ylim[0]
        else:  # "center"
            labelva = "center"
            if self.is_bigstyle:
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
        # Add offset to start & end
        feature.start += offset
        feature.end += offset
        # Add offset to exon regions
        if feature.exon_regions is not None:
            exon_regions = []
            for (exon_start, exon_end) in feature.exon_regions:
                exon_regions.append((exon_start + offset, exon_end + offset))
            feature.exon_regions = exon_regions
        return feature
