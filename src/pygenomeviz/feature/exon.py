from __future__ import annotations

from copy import deepcopy
from typing import Any

from matplotlib.figure import Axes
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from pygenomeviz.feature.base import Feature


class ExonFeature(Feature):
    """Exon Feature DataClass"""

    def __init__(
        self,
        exon_regions: list[tuple[int, int]],
        strand: int,
        label: str = "",
        labelsize: int = 15,
        labelcolor: str = "black",
        plotstyle: str = "bigarrow",  # "(big)arrow", "(big)box", "(big)rbox"
        facecolor: str = "orange",
        edgecolor: str = "black",
        linewidth: float = 0,
        labelrotation: int = 45,
        labelvpos: str = "strand",  # "top", "center", "bottom", "strand"
        labelhpos: str = "center",  # "left", "center", "right"
        labelha: str = "left",  # "left", "center", "right"
        arrow_shaft_ratio: float = 0.5,
        size_ratio: float = 0.9,
        exon_labels: list[str] | None = None,
        exon_label_kws: dict[str, Any] | None = None,
        patch_kws: dict[str, Any] | None = None,
        intron_patch_kws: dict[str, Any] | None = None,
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
