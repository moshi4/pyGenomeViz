from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from typing import Any, get_args

from Bio.SeqFeature import SeqFeature
from matplotlib.axes import Axes
from matplotlib.patches import FancyArrow, PathPatch, Rectangle
from matplotlib.path import Path

from pygenomeviz.config import LiteralTypes


@dataclass
class Feature:
    """Feature DataClass"""

    start: int
    end: int
    strand: LiteralTypes.STRAND = 1
    label: str = ""
    labelsize: int = 15
    labelcolor: str = "black"
    plotstyle: LiteralTypes.PLOTSTYLE = "bigarrow"
    facecolor: str = "orange"
    edgecolor: str = "black"
    linewidth: float = 0
    labelrotation: int = 45
    labelvpos: LiteralTypes.LABELVPOS = "strand"
    labelhpos: LiteralTypes.LABELHPOS = "center"
    labelha: LiteralTypes.LABELHA = "left"
    arrow_shaft_ratio: float = 0.5
    size_ratio: float = 0.9
    patch_kws: dict[str, Any] | None = None
    seq_feature: SeqFeature | None = None

    def __post_init__(self):
        # start-end value for HTML display
        self._display_start = self.start
        self._display_end = self.end
        # Change unknown strand value to 1
        if self.strand not in (1, -1):
            self.strand = 1
        # Check start, end postion
        if self.start > self.end:
            err_msg = f"Feature 'end' must be larger than 'start' ({self})"
            raise ValueError(err_msg)
        # Check labelvpos
        if self.labelvpos not in get_args(LiteralTypes.LABELVPOS):
            raise ValueError(f"{self.labelvpos=} is invalid parameter.")
        # Check labelhpos
        if self.labelhpos not in get_args(LiteralTypes.LABELHPOS):
            raise ValueError(f"{self.labelhpos} is invalid parameter.")
        # Check labelha
        if self.labelha not in get_args(LiteralTypes.LABELHA):
            raise ValueError(f"{self.labelha=} is invalid parameter.")
        # Check feature plot style
        plotstyle_list = get_args(LiteralTypes.PLOTSTYLE)
        if self.plotstyle not in plotstyle_list:
            err_msg = f"'plotstyle must be '{'|'.join(plotstyle_list)}'.\n"
            err_msg += f"{self.plotstyle=} is invalid."
            raise ValueError(err_msg)
        # Check arrow shaft ratio
        if not 0 <= self.arrow_shaft_ratio <= 1:
            err_msg = "'arrow_shaft_ratio' must be '0 <= value <= 1' "
            err_msg += f"({self.arrow_shaft_ratio=})"
            raise ValueError(err_msg)
        # Check size ratio
        if not 0 <= self.size_ratio <= 1:
            err_msg = f"'size_ratio' must be '0 <= value <= 1' ({self.size_ratio=})"
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
            raise ValueError(f"{self.plotstyle=} is invalid.")

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

    @property
    def gid(self) -> str:
        """Group ID"""
        start, end = self._display_start, self._display_end
        type, gene, protein_id, product = "na", "na", "na", "na"
        if self.seq_feature is not None:
            type = self.seq_feature.type
            qualifiers = self.seq_feature.qualifiers
            gene = qualifiers.get("gene", ["na"])[0]
            protein_id = qualifiers.get("protein_id", ["na"])[0].split(".")[0]
            product = qualifiers.get("product", ["na"])[0]
            if product == "na":
                product = qualifiers.get("Name", ["na"])[0]
            # Replace special characters to underscore (For html ID tag selection)
            trans_dict = {e: "_" for e in list(" /:;()+.,'`\"\\!|^~[]{}<>#$%&@?=")}
            trans_table = str.maketrans(trans_dict)
            gene = gene.translate(trans_table)
            protein_id = protein_id.translate(trans_table)
            product = product.translate(trans_table)

        return (
            f"Feature_{start}_{end}_{self.strand}_"
            + f"type_{type}_gene_{gene}_protein_id_{protein_id}_product_{product}"
        )

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
            "gid": self.gid,
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


class ExonFeature(Feature):
    """Exon Feature DataClass"""

    def __init__(
        self,
        exon_regions: list[tuple[int, int]],
        strand: LiteralTypes.STRAND = 1,
        label: str = "",
        labelsize: int = 15,
        labelcolor: str = "black",
        plotstyle: LiteralTypes.PLOTSTYLE = "bigarrow",
        facecolor: str = "orange",
        edgecolor: str = "black",
        linewidth: float = 0,
        labelrotation: int = 45,
        labelvpos: LiteralTypes.LABELVPOS = "strand",
        labelhpos: LiteralTypes.LABELHPOS = "center",
        labelha: LiteralTypes.LABELHA = "left",
        arrow_shaft_ratio: float = 0.5,
        size_ratio: float = 0.9,
        exon_labels: list[str] | None = None,
        exon_label_kws: dict[str, Any] | None = None,
        patch_kws: dict[str, Any] | None = None,
        intron_patch_kws: dict[str, Any] | None = None,
        seq_feature: SeqFeature | None = None,
    ):
        self.exon_regions = exon_regions
        self.exon_labels = exon_labels
        self.exon_label_kws = {} if exon_label_kws is None else exon_label_kws
        self.intron_patch_kws = {} if intron_patch_kws is None else intron_patch_kws
        self._check_exon_regions()
        super().__init__(
            exon_regions[0][0],  # start
            exon_regions[-1][-1],  # end
            strand,
            label,
            labelsize,
            labelcolor,
            plotstyle,
            facecolor,
            edgecolor,
            linewidth,
            labelrotation,
            labelvpos,
            labelhpos,
            labelha,
            arrow_shaft_ratio,
            size_ratio,
            patch_kws,
            seq_feature,
        )

        if exon_labels is not None and len(exon_regions) != len(exon_labels):
            err_msg = "'exon_regions' & 'exon_labels' length is diffrent."
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

        # Plot intron feature
        for (start, end) in self.intron_regions:
            p = self._intron_patch(start, end, ylim)
            ax.add_patch(p)

        # Plot exon feature
        exon_regions = self.exon_regions[:: self.strand]
        for idx, (start, end) in enumerate(exon_regions, 1):
            if self.plotstyle in ("bigbox", "box"):
                p = self._box_patch(start, end, ylim)
            elif self.plotstyle in ("bigarrow", "arrow"):
                no_head_length = False if idx == len(exon_regions) else True
                p = self._arrow_patch(start, end, ylim, max_track_size, no_head_length)
            elif self.plotstyle in ("bigrbox", "rbox"):
                p = self._rbox_patch(start, end, ylim, max_track_size)
            else:
                raise ValueError(f"{self.plotstyle=} is invalid.")
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
        super().plot_label(ax, ylim)
        if self.exon_labels is not None:
            for (start, end), label in zip(self.exon_regions, self.exon_labels):
                if label != "" and self.labelsize != 0:
                    label_kwargs = self._label_kwargs(start, end, label, ylim)
                    label_kwargs.update(self.exon_label_kws)
                    ax.text(**label_kwargs)

    @property
    def intron_regions(self) -> list[tuple[int, int]]:
        """Intron regions"""
        intron_regions = []
        if len(self.exon_regions) > 1:
            for i in range(0, len(self.exon_regions) - 1):
                intron_start = self.exon_regions[i][1]
                intron_end = self.exon_regions[i + 1][0]
                if intron_end < intron_start:
                    continue
                intron_regions.append((intron_start, intron_end))
        return intron_regions

    def _check_exon_regions(self) -> None:
        """Check exon_regions values are properly set"""
        max_pos_record = None
        for (exon_start, exon_end) in self.exon_regions:
            if exon_start >= exon_end:
                err_msg = "Exon start-end value must be 'start < end'."
                raise ValueError(err_msg)
            if max_pos_record is not None and exon_start < max_pos_record:
                err_msg = "Set exon regions in ascending order of genomic position."
                raise ValueError(err_msg)
            max_pos_record = exon_end

    def _intron_patch(
        self, start: int, end: int, ylim: tuple[float, float]
    ) -> PathPatch:
        """Intron patch

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
        intron_patch : PathPatch
            Intron patch
        """
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
        return PathPatch(Path(verts, codes), **self._intron_patch_kwargs())

    def _intron_patch_kwargs(self) -> dict[str, Any]:
        """Intron patch keyword arguments dict

        Returns
        -------
        intron_patch_kwargs : dict[str, Any]
            Intron patch keyword arguments dict
        """
        return {
            "lw": 1,
            "fill": False,
            "clip_on": False,
            "zorder": 5 if self.is_bigstyle else -5,
            **self.intron_patch_kws,
        }

    def __add__(self, offset: int):
        feature = deepcopy(self)
        # Add offset to start & end
        feature.start += offset
        feature.end += offset
        # Add offset to exon regions
        exon_regions = []
        for (exon_start, exon_end) in feature.exon_regions:
            exon_regions.append((exon_start + offset, exon_end + offset))
        feature.exon_regions = exon_regions
        return feature
