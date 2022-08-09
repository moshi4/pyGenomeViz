from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from typing import Any

from matplotlib.figure import Axes
from matplotlib.patches import FancyArrow, PathPatch, Rectangle
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
    patch_kws: dict[str, Any] | None = None

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

    def plot_feature(
        self, ax: Axes, max_track_size: int, ylim: tuple[float, float]
    ) -> None:
        """Plot feature

        Parameters
        ----------
        ax : Axes
            Matplotlib axes object to be plotted
        max_track_size : int
            Max track size
        ylim : tuple[float, float]
            Y-axis limit
        """
        ylim = (ylim[0] * self.size_ratio, ylim[1] * self.size_ratio)
        start, end = self.start, self.end

        if self.plotstyle in ("bigbox", "box"):
            p = self._box_patch(start, end, ylim)
        elif self.plotstyle in ("bigarrow", "arrow"):
            p = self._arrow_patch(start, end, ylim, max_track_size)
        elif self.plotstyle in ("bigrbox", "rbox"):
            p = self._rbox_patch(start, end, ylim, max_track_size)
        else:
            raise ValueError(f"'{self.plotstyle}' is invalid plotstyle.")

        ax.add_patch(p)

    def plot_label(self, ax: Axes, ylim: tuple[float, float]) -> None:
        """Plot label

        Parameters
        ----------
        ax : Axes
            Matplotlib axes object to be plotted
        ylim : tuple[float, float]
            Y-axis limit
        """
        if self.label != "" and self.labelsize != 0:
            start, end = self.start, self.end
            ax.text(**self._label_kwargs(start, end, self.label, ylim))

    @property
    def length(self) -> int:
        """Feature length"""
        return self.end - self.start

    @property
    def is_bigstyle(self) -> bool:
        """Check plotstyle is 'big~~~' or not"""
        return self.plotstyle.startswith("big")

    def _box_patch(self, start: int, end: int, ylim: tuple[float, float]) -> Rectangle:
        """Box patch

        Parameters
        ----------
        start : int
            Start position
        end : int
            End position
        ylim : tuple[float, float]
            Y-axis limit

        Returns
        -------
        box_patch : Rectangle
            Box patch
        """
        # x, y
        x = start
        if self.is_bigstyle or self.strand == -1:
            y = ylim[0]
        else:
            y = 0
        # width, height
        width = end - start
        if self.is_bigstyle:
            height = ylim[1] - ylim[0]
        else:
            height = ylim[1]

        return Rectangle((x, y), width, height, **self._patch_kwargs())

    def _arrow_patch(
        self,
        start: int,
        end: int,
        ylim: tuple[float, float],
        max_track_size: int,
        no_head_length: bool = False,
    ) -> FancyArrow:
        """Arrow patch

        Parameters
        ----------
        start : int
            Start position
        end : int
            End position
        ylim : tuple[float, float]
            Y-axis limit
        max_track_size : int
            Max track size (Use for head length calculation)
        no_head_length : bool, optional
            If True, set head length as 0

        Returns
        -------
        arrow_patch : FancyArrow
            Arrow patch
        """
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
        length = end - start
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
        if length < head_length:
            head_length = length

        if no_head_length:
            head_length = 0
            head_width = shaft_width

        return FancyArrow(
            x,
            y,
            dx,
            dy,
            width=shaft_width,
            length_includes_head=True,
            head_width=head_width,
            head_length=head_length,
            **self._patch_kwargs(),
        )

    def _rbox_patch(
        self, start: int, end: int, ylim: tuple[float, float], max_track_size: int
    ) -> PathPatch:
        """Rounded box patch

        Parameters
        ----------
        start : int
            Start position
        end : int
            End position
        ylim : tuple[float, float]
            Y-axis limit
        max_track_size : int
            Max track size (Use for rounded size calculation)

        Returns
        -------
        rbox_patch : PathPatch
            Rounded box patch
        """
        length = end - start
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
        return PathPatch(Path(verts, codes), **self._patch_kwargs())

    def _patch_kwargs(self) -> dict[str, Any]:
        """Patch keyword arguments dict

        Returns
        -------
        patch_kwargs : dict[str, Any]
            Patch keyword arguments dict
        """
        patch_kws = {} if self.patch_kws is None else self.patch_kws
        return {
            "fc": self.facecolor,
            "ec": self.edgecolor,
            "lw": self.linewidth,
            "clip_on": False,
            "zorder": 5 if self.is_bigstyle else -5,
            **patch_kws,
        }

    def _label_kwargs(
        self, start: int, end: int, label: str, ylim: tuple[float, float]
    ) -> dict[str, Any]:
        """Label keyword arguments dict

        Parameters
        ----------
        start : int
            Start position
        end : int
            End position
        label : str
            Label
        ylim : tuple[float, float]
            Y-axis limit

        Returns
        -------
        label_kwargs : dict[str, Any]
            Label keyword arguments dict
        """
        # x
        if self.labelhpos == "left":
            x = start
        elif self.labelhpos == "right":
            x = end
        else:  # "center"
            x = (start + end) / 2

        # labelrotation
        labelrotation = self.labelrotation
        if self.labelvpos == "strand":
            labelrotation = labelrotation * self.strand

        # labelvpos
        labelvpos = self.labelvpos
        if labelvpos == "strand":
            labelvpos = "bottom" if self.strand == -1 else "top"

        # labelva, y
        label_margin = 0.1
        ylim = (ylim[0] * self.size_ratio, ylim[1] * self.size_ratio)
        if labelvpos == "top":
            labelva = "bottom"
            y = ylim[1] + label_margin
        elif labelvpos == "bottom":
            labelva = "top"
            y = ylim[0] - label_margin
        else:  # "center"
            labelva = "center"
            if self.is_bigstyle:
                y = 0
            else:
                y = (abs(ylim[0]) / 2) * self.strand

        # labelrotation=90, labelha="center", rotation_mode="default" is best
        rotation_mode = "default" if self.labelrotation == 90 else "anchor"

        return {
            "x": x,
            "y": y,
            "s": label,
            "color": self.labelcolor,
            "fontsize": self.labelsize,
            "rotation": labelrotation,
            "ha": self.labelha,
            "va": labelva,
            "zorder": 10,
            "rotation_mode": rotation_mode,
        }

    def __add__(self, offset: int):
        feature = deepcopy(self)
        feature.start += offset
        feature.end += offset
        return feature
