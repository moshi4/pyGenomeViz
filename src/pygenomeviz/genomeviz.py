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
        tick_track_ratio: float = 1.0,
        track_spines: bool = False,
        tick_style: Optional[str] = None,
        plot_size_thr: float = 0.0005,  # 0.05 %
        tick_labelsize: int = 15,
    ):
        """
        Parameters
        ----------
        fig_width : float, optional
            Figure width
        fig_track_height : float, optional
            Figure track height
        align_type : str, optional
            Track align type
        feature_size_ratio : float, optional
            Feature size ratio to track
        link_size_ratio : float, optional
            Link size ratio to track
        arrow_shaft_ratio : float, optional
            Feature arrow shaft ratio
        feature_track_ratio : float, optional
            Feature track ratio
        link_track_ratio : float, optional
            Link track ratio
        tick_track_ratio : float, optional
            Tick track ratio
        track_spines : bool, optional
            Display track spines
        tick_style : Optional[str], optional
            Tick style (`axis`|`bar`)
        plot_size_thr : float, optional
            Plot size threshold
        tick_labelsize : int, optional
            Tick label size

        Notes
        -----
        If `plot_size_thr=0.0005` and `max_track_size=4.0Mb`, features smaller than
        `max_track_size * plot_size_thr=2.0Kb` are not plotted.
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
        self.tick_track_ratio = tick_track_ratio
        self.tick_style = tick_style
        self.plot_size_thr = plot_size_thr
        self.tick_labelsize = tick_labelsize
        self._tracks: List[Track] = []

        self._check_init_values()

    def _check_init_values(self) -> None:
        """Check initialize values"""
        if self.align_type not in ("left", "center", "right"):
            err_msg = f"Invalid align type '{self.align_type}'."
            raise ValueError(err_msg)

        if self.tick_style is not None and self.tick_style not in ("axis", "bar"):
            err_msg = f"Invalid tick type '{self.tick_style}'."
            raise ValueError(err_msg)

        range_check_dict = {
            "feature_size_ratio": self.feature_size_ratio,
            "link_size_ratio": self.link_size_ratio,
            "arrow_shaft_ratio": self.arrow_shaft_ratio,
        }
        err_msg = ""
        for k, v in range_check_dict.items():
            if not 0 <= v <= 1:
                err_msg += f"'{k}' must be '0 <= value <= 1' (value={v})\n"
        if err_msg:
            raise ValueError(err_msg)

    @property
    def track_num(self) -> int:
        """Track number"""
        return len(self._tracks)

    @property
    def max_track_size(self) -> int:
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

    def _get_track_idx(self, track_name: str) -> int:
        """Get track index from track name

        Parameters
        ----------
        track_name : str
            Track name

        Returns
        -------
        track_idx : int
            Track index
        """
        return [t.name for t in self._tracks].index(track_name)

    def _get_link_track(
        self, feature_track_name1: str, feature_track_name2: str
    ) -> LinkTrack:
        """Get link track from two adjacent feature track names

        Parameters
        ----------
        feature_track_name1 : str
            Feature track name1
        feature_track_name2 : str
            Feature track name2

        Returns
        -------
        link_track : LinkTrack
            Target link track
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

        Parameters
        ----------
        track : Track
            Target track offset for alignment

        Returns
        -------
        track_offset : int
            Track offset
        """
        if self.align_type == "left":
            return 0
        elif self.align_type == "center":
            return int((self.max_track_size - track.size) / 2)
        elif self.align_type == "right":
            return self.max_track_size - track.size
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

        Parameters
        ----------
        name : str
            Track name
        size : int
            Track size
        labelsize : int, optional
            Track label size
        linewidth : int, optional
            Trakc line width

        Returns
        -------
        feature_track : FeatureTrack
            Feature track
        """
        # Check specified track name is unique or not
        if name in [t.name for t in self._tracks]:
            err_msg = f"track.name='{name}' is already exists. Change to another name."
            raise ValueError(err_msg)
        # Add link track between feature tracks
        if len(self._tracks) > 0:
            link_track = LinkTrack(
                f"{self._tracks[-1].name}-{name}",
                self.track_spines,
                self.link_track_ratio,
            )
            self._tracks.append(link_track)
        # Add feature track
        feature_track = FeatureTrack(
            name,
            size,
            labelsize,
            linewidth,
            self.track_spines,
            self.feature_track_ratio,
        )
        self._tracks.append(feature_track)
        return feature_track

    def add_link(
        self,
        track_link1: Tuple[str, int, int],
        track_link2: Tuple[str, int, int],
        normal_color: str = "grey",
        inverted_color: str = "red",
        interpolation_value: Optional[float] = None,
        vmin: float = 0,
        vmax: float = 100,
    ) -> None:
        """Add link data to link track

        Parameters
        ----------
        track_link1 : Tuple[str, int, int]
            Track link1 (track_name, start, end)
        track_link2 : Tuple[str, int, int]
            Track link2 (track_name, start, end)
        normal_color : str, optional
            Normal link color
        inverted_color : str, optional
            Inverted link color
        interpolation_value : Optional[float], optional
            Value for color interpolation
        vmin : float, optional
            Min value for color interpolation
        vmax : float, optional
            Max value for color interpolation
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
                interpolation_value,
                vmin,
                vmax,
            )
        )

    def plotfig(self) -> Figure:
        """Plot figure

        Returns
        -------
        figure : Figure
            Plot figure result
        """
        if self.track_num == 0:
            raise ValueError("No tracks are defined for plotting figure.")
        if self.tick_style is not None:
            self._tracks.append(
                TickTrack(
                    self.max_track_size,
                    self.tick_labelsize,
                    self.track_spines,
                    self.tick_track_ratio,
                    self.tick_style,
                )
            )

        figsize = (self.fig_width, self.fig_track_height * self.track_num)
        figure: Figure = plt.figure(figsize=figsize, facecolor="white")
        track_ratios = [t.ratio for t in self._tracks]
        spec = gridspec.GridSpec(
            nrows=self.track_num, ncols=1, height_ratios=track_ratios, hspace=0
        )
        for idx, track in enumerate(self._tracks):
            # Create new track subplot
            xlim, ylim = (0, self.max_track_size), track.ylim
            ax: Axes = figure.add_subplot(
                spec[idx], xlim=xlim, ylim=ylim, fc="none", zorder=track.zorder
            )
            # Set 'spines' and 'ticks' visibility
            for spine, display in track.spines_params.items():
                ax.spines[spine].set_visible(display)
            ax.tick_params(**track.tick_params)

            if isinstance(track, FeatureTrack):
                # Plot track label
                if track.labelsize != 0:
                    ax.text(-self.max_track_size * 0.01, 0, **track.label_params)
                # Plot track scale line
                track_offset = self._track_offset(track)
                xmin, xmax = track_offset, track.size + track_offset
                ax.hlines(0, xmin, xmax, "black", linewidth=track.linewidth)

                for feature in [f + track_offset for f in track.features]:
                    # Don't draw too small feature (To reduce drawing time)
                    if feature.length < self.max_track_size * self.plot_size_thr:
                        continue
                    # Plot feature object
                    obj_params = feature.obj_params(
                        self.max_track_size,
                        ylim,
                        self.feature_size_ratio,
                        self.arrow_shaft_ratio,
                    )
                    ax.arrow(**obj_params)
                    # Plot feature text
                    if feature.labelsize != 0:
                        ax.text(**feature.text_params(ylim, self.feature_size_ratio))

            elif isinstance(track, LinkTrack):
                for link in track.links:
                    link = link.add_offset(self._track_name2offset)
                    link_ymin = ylim[0] * self.link_size_ratio
                    link_ymax = ylim[1] * self.link_size_ratio
                    xy = link.polygon_xy(link_ymin, link_ymax)
                    p = patches.Polygon(xy=xy, fc=link.color, ec=link.color)
                    ax.add_patch(p)

            elif isinstance(track, TickTrack):
                if self.tick_style == "axis":
                    ax.xaxis.set_major_locator(MaxNLocator(10, steps=[1, 2, 5, 10]))
                    ax.xaxis.set_major_formatter(track.tick_formatter)
                elif self.tick_style == "bar":
                    ymin, ycenter, ymax = track.ymin, track.ycenter, track.ymax
                    ax.hlines(ycenter, track.xmin, track.xmax, "black", linewidth=1)
                    ax.vlines(track.xmin, ymin, ymax, "black", linewidth=1)
                    ax.vlines(track.xmax, ymin, ymax, "black", linewidth=1)
                    ax.text(**track.scalebar_text_params)

        return figure

    def savefig(
        self,
        savefile: Union[str, Path, BytesIO],
        dpi: int = 300,
        pad_inches: float = 0.5,
    ) -> None:
        """Save figure to file

        Parameters
        ----------
        savefile : Union[str, Path, BytesIO]
            Save file
        dpi : int, optional
            DPI
        pad_inches : float, optional
            Padding inches
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
