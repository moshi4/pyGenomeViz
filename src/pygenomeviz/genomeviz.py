from __future__ import annotations

from io import BytesIO
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
from matplotlib import gridspec, patches
from matplotlib.figure import Axes, Figure
from matplotlib.ticker import MaxNLocator

from pygenomeviz.link import Link
from pygenomeviz.track import FeatureTrack, LinkTrack, TickTrack, Track


class GenomeViz:
    """GenomeViz Class"""

    def __init__(
        self,
        fig_width: float = 15,
        fig_track_height: float = 1.0,
        align_type: str = "left",
        feature_size_ratio: float = 0.9,
        link_size_ratio: float = 0.9,
        arrow_shaft_ratio: float = 0.5,
        feature_track_ratio: float = 1.0,
        link_track_ratio: float = 1.0,
        track_spines: bool = False,
        tick_type: Optional[str] = None,
    ):
        """GenomeViz constructor

        Args:
            fig_width (float, optional): Figure width
            fig_track_height (float, optional): Figure track height
            align_type (str, optional): Track align type ('left'|'center'|'right')
            feature_size_ratio (float, optional): Feature size ratio  [0.0 - 1.0]
            link_size_ratio (float, optional): Link size ratio [0.0 - 1.0]
            arrow_shaft_ratio (float, optional): Feature arrow shaft ratio [0.0 - 1.0]
            feature_track_ratio (float, optional): Feature track ratio
            link_track_ratio (float, optional): Link track ratio
            track_spines (bool, optional): Display track spines
            tick_type (Optional[str]): Tick type ('all' or 'partial')
        """
        self.fig_width = fig_width
        self.fig_track_height = fig_track_height
        self.align_type = align_type
        self.track_spines = track_spines
        self.feature_size_ratio = feature_size_ratio
        self.link_size_ratio = link_size_ratio
        self.arrow_shaft_ratio = arrow_shaft_ratio
        self.feature_track_ratio = feature_track_ratio
        self.link_track_ratio = link_track_ratio
        self.tick_type = tick_type
        self._min_plot_size = 0.001  # 0.1 %
        self._tracks: List[Track] = []
        self._tick_labelsize = 15

        self._check_arg_values()

    def _check_arg_values(self) -> None:
        if self.align_type not in ("left", "center", "right"):
            err_msg = f"Invalid align type '{self.align_type}'."
            raise ValueError(err_msg)

        if self.tick_type is not None and self.tick_type not in ("all", "partial"):
            err_msg = f"Invalid tick type '{self.tick_type}'."
            raise ValueError(err_msg)

        range_check_dict = {
            "feature_size_ratio": self.feature_size_ratio,
            "link_size_ratio": self.link_size_ratio,
            "arrow_shaft_ratio": self.arrow_shaft_ratio,
            "feature_track_ratio": self.feature_track_ratio,
            "link_track_ratio": self.link_track_ratio,
        }
        err_msg = ""
        for k, v in range_check_dict.items():
            if not 0 <= v <= 1:
                err_msg += f"'{k}' must be '0 <= value <= 1' ({k}={v})\n"
        if err_msg:
            raise ValueError(err_msg)

    @property
    def _track_num(self) -> int:
        """Track number"""
        return len(self._tracks)

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
                track_ratios.append(self.feature_track_ratio)
            elif isinstance(track, LinkTrack):
                track_ratios.append(self.link_track_ratio)
            elif isinstance(track, TickTrack):
                track_ratios.append(self.feature_track_ratio * 0.5)
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
        """Get link track by two adjacent feature track names

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
            else:
                err_msg = "LinkTrack is not found between target feature tracks "
                err_msg += f"(`{feature_track_name1}` and '{feature_track_name2}')"
                raise ValueError(err_msg)
        else:
            err_msg = f"`{feature_track_name1}` and '{feature_track_name2}' "
            err_msg += "are not adjacent feature tracks."
            raise ValueError(err_msg)

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
        labelsize: int = 20,
        linewidth: int = 1,
    ) -> FeatureTrack:
        """Add feature track

        Args:
            name (str): Track name
            size (int): Track size
            labelsize (int, optional): Track label size
            linewidth (int, optional): Track line width

        Returns:
            FeatureTrack: Newly added FeatureTrack
        """
        # Check specified track name is unique or not
        if name in [t.name for t in self._tracks]:
            err_msg = f"track.name='{name}' is already exists. Change to another name."
            raise ValueError(err_msg)
        # Add link track between feature tracks
        if len(self._tracks) > 0:
            link_track = LinkTrack(f"{self._tracks[-1].name}-{name}", self.track_spines)
            self._tracks.append(link_track)
        # Add feature track
        feature_track = FeatureTrack(
            name, size, labelsize, linewidth, self.track_spines
        )
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
            track_link1 (Tuple[str, int, int]): Track link1 (track name, start, end)
            track_link2 (Tuple[str, int, int]): Track link2 (track name, start, end)
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
                above_track_link[0],
                above_track_link[1],
                above_track_link[2],
                below_track_link[0],
                below_track_link[1],
                below_track_link[2],
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
        if self._track_num == 0:
            raise ValueError("No tracks are defined for plotting figure.")
        if self.tick_type is not None:
            self._tracks.append(
                TickTrack(
                    self._max_track_size,
                    self._tick_labelsize,
                    self.track_spines,
                    self.tick_type,
                )
            )

        figsize = (self.fig_width, self.fig_track_height * self._track_num)
        figure: Figure = plt.figure(figsize=figsize, tight_layout=False)
        spec = gridspec.GridSpec(
            nrows=self._track_num, ncols=1, height_ratios=self._track_ratios, hspace=0
        )
        for idx, track in enumerate(self._tracks):
            # Create new track subplot
            xlim, ylim = (0, self._max_track_size), (-1.0, 1.0)
            ax: Axes = figure.add_subplot(spec[idx], xlim=xlim, ylim=ylim)

            # Set 'spines' and 'ticks' visibility
            for spine, display in track.spines_params.items():
                ax.spines[spine].set_visible(display)
            ax.tick_params(**track.tick_params)

            if isinstance(track, FeatureTrack):
                # Plot track label
                if track.labelsize != 0:
                    ax.text(-self._max_track_size * 0.01, 0, **track.label_params)
                # Plot track scale line
                track_offset = self._track_offset(track)
                xmin, xmax = track_offset, track.size + track_offset
                ax.hlines(0, xmin, xmax, "black", linewidth=track.linewidth)

                for feature in [f + track_offset for f in track.features]:
                    # Plot feature object
                    obj_params = feature.obj_params(
                        self._max_track_size,
                        ylim,
                        self.feature_size_ratio,
                        self.arrow_shaft_ratio,
                    )
                    ax.arrow(**obj_params)
                    # Plot feature text
                    feature_text_x = (feature.start + feature.end) / 2
                    if feature.labelsize != 0:
                        ax.text(feature_text_x, obj_params["y"], **feature.text_params)

            elif isinstance(track, LinkTrack):
                for link in track.links:
                    link = link.add_offset(self._track_name2offset)
                    link_ymin = ylim[0] * self.link_size_ratio
                    link_ymax = ylim[1] * self.link_size_ratio
                    xy = link.polygon_xy(link_ymin, link_ymax)
                    p = patches.Polygon(xy=xy, fc=link.color, ec=link.color)
                    ax.add_patch(p)

            elif isinstance(track, TickTrack):
                if self.tick_type == "all":
                    ax.xaxis.set_major_locator(MaxNLocator(10, steps=[1, 2, 5, 10]))
                    ax.xaxis.set_major_formatter(track.tick_formatter)
                elif self.tick_type == "partial":
                    ax.hlines(ylim[0] / 2, track.xmin, track.xmax, "black", linewidth=1)
                    ax.vlines(track.xmin, ylim[0], 0, "black", linewidth=1)
                    ax.vlines(track.xmax, ylim[0], 0, "black", linewidth=1)
                    x, y, label = track.xcenter, ylim[0], track.scalebar_label
                    ax.text(x, y, label, fontsize=15, ha="center", va="top")

        return figure

    def savefig(
        self,
        savefile: Union[str, Path, BytesIO],
        dpi: int = 300,
        pad_inches: float = 0.5,
    ) -> None:
        """Save figure to file

        Args:
            savefile (Union[str, Path, BytesIO]): Save file
            dpi (int, optional): DPI
            pad_inches (float, optional): Padding inches
        """
        figure = self.plotfig()
        figure.savefig(
            fname=savefile,
            dpi=dpi,
            pad_inches=pad_inches,
            bbox_inches="tight",
        )

    def print_tracks_info(self) -> None:
        """Print tracks info (For developer debugging)"""
        for idx, track in enumerate(self._tracks, 1):
            class_name = track.__class__.__name__
            print(f"\n# Track{idx:02d}: Name='{track.name}' ({class_name})")
            if isinstance(track, FeatureTrack):
                print(f"# Size={track.size}, FeatureNumber={len(track.features)}")
                for feature in track.features:
                    print(feature)
            elif isinstance(track, LinkTrack):
                print(f"# Size={track.size}, LinkNumber={len(track.links)}")
                for link in track.links:
                    print(link)
            elif isinstance(track, TickTrack):
                pass
