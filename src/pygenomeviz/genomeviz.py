from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
from matplotlib import gridspec, patches
from matplotlib.figure import Axes, Figure

from pygenomeviz.link import Link
from pygenomeviz.track import FeatureTrack, LinkTrack, Track


class GenomeViz:
    def __init__(
        self,
        fig_width: float = 15,
        fig_track_height: float = 1.0,
        align_type: str = "left",
    ):
        """GenomeViz constructor

        Args:
            fig_width (float, optional): Figure width. Defaults to 15.
            fig_track_height (float, optional): Figure track height. Defaults to 1.0.
            align_type (str, optional): Track align type. Defaults to "left".
        """
        self.fig_width = fig_width
        self.fig_track_height = fig_track_height
        self.align_type = align_type
        self._min_plot_size = 0.001  # 0.1 %
        self._feature_track_ratio = 1.0
        self._link_track_ratio = 1.0
        self._track_y_pad = 0.05
        self._feature_size_ratio = 1.0
        self._tracks: List[Union[FeatureTrack, LinkTrack]] = []

        if self.align_type not in ("left", "center", "right"):
            err_msg = f"Invalid align type '{self.align_type}'."
            raise ValueError(err_msg)

    @property
    def _max_track_size(self) -> int:
        """Max track size"""
        if len(self._tracks) > 0:
            return max([track.size for track in self._tracks])
        else:
            return 0

    @property
    def _track_name2offset(self) -> Dict[str, int]:
        """Track name & offset dict"""
        track_name2offset = {}
        for track in self._tracks:
            track_name2offset[track.name] = self._track_offset(track)
        return track_name2offset

    @property
    def _track_ratios(self) -> List[float]:
        """Each track height ratio list"""
        track_ratios = []
        for track in self._tracks:
            if isinstance(track, FeatureTrack):
                track_ratios.append(self._feature_track_ratio)
            elif isinstance(track, LinkTrack):
                track_ratios.append(self._link_track_ratio)
        return track_ratios

    def _get_track_idx(self, track_name: str) -> int:
        """Get track index by track name

        Args:
            track_name (str): Track name

        Returns:
            int: Track index
        """
        track_names = [t.name for t in self._tracks]
        return track_names.index(track_name)

    def _get_link_track(
        self, feature_track_name1: str, feature_track_name2: str
    ) -> LinkTrack:
        """Get link track by two feature track names

        Args:
            feature_track_name1 (str): Feature track name1
            feature_track_name2 (str): Feature track name2

        Returns:
            LinkTrack: Target LinkTrack
        """
        feature_track_idx1 = self._get_track_idx(feature_track_name1)
        feature_track_idx2 = self._get_track_idx(feature_track_name2)
        if abs(feature_track_idx1 - feature_track_idx2) == 2:
            link_track_idx = int((feature_track_idx1 + feature_track_idx2) / 2)
            target_link_track = self._tracks[link_track_idx]
            if isinstance(target_link_track, LinkTrack):
                return target_link_track

        # TODO: Error handling
        raise RuntimeError()

    def _track_offset(self, track: Track) -> int:
        """Get track offset

        Args:
            track (Track): Target track for alignment offset calculation

        Returns:
            int: Target track offset
        """
        if self.align_type == "left":
            return 0
        elif self.align_type == "center":
            return int((self._max_track_size - track.size) / 2)
        elif self.align_type == "right":
            return self._max_track_size - track.size
        else:
            return 0

    def add_feature_track(
        self,
        name: str,
        size: int,
        labelsize: int = 30,
        linewidth: int = 2,
    ) -> FeatureTrack:
        """Add feature track

        Args:
            name (str): Track name
            size (int): Track size
            labelsize (int, optional): Track label size. Defaults to 30.
            linewidth (int, optional): Track line width. Defaults to 2.

        Returns:
            FeatureTrack: Newly added FeatureTrack
        """
        if len(self._tracks) > 0:
            link_track = LinkTrack(f"{self._tracks[-1].name}-{name}", 0)
            self._tracks.append(link_track)
        feature_track = FeatureTrack(name, size, labelsize, linewidth)
        self._tracks.append(feature_track)
        return feature_track

    def add_link(
        self,
        track_link1: Tuple[str, int, int],
        track_link2: Tuple[str, int, int],
        normal_color: str = "grey",
        inverted_color: str = "red",
        identity: Optional[float] = None,
        interpolation: bool = True,
    ) -> None:
        """Add link data to link track

        Args:
            track_link1 (Tuple[str, int, int]): Track link data1
            track_link2 (Tuple[str, int, int]): Track link data2
            normal_color (str, optional): Normal link color.
            inverted_color (str, optional): Inverted link color.
            identity (Optional[float], optional): Link identity [0 - 100].
            interpolation (bool, optional): Enable color interpolation by identity.
        """
        link_track = self._get_link_track(track_link1[0], track_link2[0])
        if self._get_track_idx(track_link1[0]) < self._get_track_idx(track_link2[0]):
            above_track_link, below_track_link = track_link1, track_link2
        else:
            above_track_link, below_track_link = track_link2, track_link2
        link_track.add_link(
            Link(
                *above_track_link,
                *below_track_link,
                normal_color,
                inverted_color,
                identity,
                interpolation,
            )
        )

    def plotfig(self) -> Figure:
        """Plot figure

        Returns:
            Figure: Plot figure result
        """
        track_num = len(self._tracks)
        if track_num == 0:
            raise RuntimeError("No tracks are defined for plotting figure.")
        figsize = (self.fig_width, self.fig_track_height * track_num)
        figure: Figure = plt.figure(figsize=figsize, tight_layout=False)
        spec = gridspec.GridSpec(
            nrows=track_num, ncols=1, height_ratios=self._track_ratios, hspace=0
        )
        for idx, track in enumerate(self._tracks):
            # Create new track subplot
            xlim = (0, self._max_track_size)
            ylim = (0 - self._track_y_pad, 1.0 + self._track_y_pad)
            y_center = (ylim[1] + ylim[0]) / 2
            ax: Axes = figure.add_subplot(spec[idx], xlim=xlim, ylim=ylim)

            # Disable 'spines' and 'ticks' visibility
            for spine in ax.spines:
                ax.spines[spine].set_visible(False)
            ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

            # Plot track label
            x = -self._max_track_size / 100
            ax.text(
                x=x,
                y=y_center,
                s=track.name,
                fontsize=track.labelsize,
                ha="right",
                va="center",
            )

            # Plot track scale line
            track_offset = self._track_offset(track)
            xmin, xmax = track_offset, track.size + track_offset
            ax.hlines(
                y_center, xmin, xmax, "black", linewidth=track.linewidth, zorder=-1
            )

            if isinstance(track, FeatureTrack):
                offset_features = [f + track_offset for f in track.features]
                for feature in offset_features:
                    x = feature.start if feature.strand == 1 else feature.end
                    arrow_length = feature.length * feature.strand

                    head_length = self._max_track_size * 0.02
                    if abs(feature.length) < head_length:
                        head_length = abs(feature.length)

                    zorder = -5
                    if feature.plotstyle == "bigarrow":
                        y, width, head_width = 0.5, 0.4, 1.0
                        zorder = 5
                    elif feature.plotstyle == "arrow":
                        y = 0.75 if feature.strand == 1 else 0.25
                        width, head_width = 0.2, 0.5
                    elif feature.plotstyle == "box":
                        y = 0.75 if feature.strand == 1 else 0.25
                        width, head_width, head_length = 0.5, 0.5, 0
                    else:
                        raise ValueError()

                    ax.arrow(
                        x=x,
                        y=y,
                        dx=arrow_length,
                        dy=0,
                        shape="full",
                        fc=feature.facecolor,
                        ec=feature.edgecolor,
                        alpha=1.0,
                        width=width * self._feature_size_ratio,
                        head_width=head_width * self._feature_size_ratio,
                        head_length=head_length,
                        length_includes_head=True,
                        zorder=zorder,
                    )

            elif isinstance(track, LinkTrack):
                for link in track.links:
                    link = link.add_offset(self._track_name2offset)
                    polygon_xy = [
                        (link.track_start2, 0),
                        (link.track_end2, 0),
                        (link.track_end1, 1),
                        (link.track_start1, 1),
                    ]
                    p = patches.Polygon(xy=polygon_xy, fc=link.color, ec=link.color)
                    ax.add_patch(p)

        return figure

    def savefig(
        self, savefile: Union[str, Path], dpi: int = 300, pad_inches: float = 0.5
    ) -> None:
        """Save figure to file

        Args:
            savefile (Union[str, Path]): Save file
            dpi (int, optional): DPI. Defaults to 300.
            pad_inches (float, optional): Padding inches. Defaults to 0.5.
        """
        figure = self.plotfig()
        figure.savefig(
            fname=savefile,
            dpi=dpi,
            pad_inches=pad_inches,
            bbox_inches="tight",
        )
