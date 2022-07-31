from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
from matplotlib import colors, gridspec
from matplotlib.colorbar import ColorbarBase
from matplotlib.figure import Axes, Figure
from matplotlib.ticker import MaxNLocator

from pygenomeviz.link import Link
from pygenomeviz.track import FeatureSubTrack, FeatureTrack, LinkTrack, TickTrack, Track

warnings.filterwarnings("ignore")


class GenomeViz:
    """GenomeViz Class"""

    def __init__(
        self,
        fig_width: float = 15,
        fig_track_height: float = 1.0,
        align_type: str = "left",
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
            (Figure height = `track number` * `fig track height`)
        align_type : str, optional
            Track align type (`left`|`center`|`height`)
        feature_track_ratio : float, optional
            Feature track ratio
        link_track_ratio : float, optional
            Link track ratio
        tick_track_ratio : float, optional
            Tick track ratio
        track_spines : bool, optional
            If True, display track spines
        tick_style : Optional[str], optional
            Tick style (`axis`|`bar`)
        plot_size_thr : float, optional
            Plot feature, link size threshold.
            If `plot_size_thr=0.0005` and `max_track_size=4.0Mb`, feature, link
            smaller than `max_track_size * plot_size_thr=2.0Kb` are not plotted.
        tick_labelsize : int, optional
            Tick label size
        """
        self.fig_width = fig_width
        self.fig_track_height = fig_track_height
        self.align_type = align_type
        self.track_spines = track_spines
        self.feature_track_ratio = feature_track_ratio
        self.link_track_ratio = link_track_ratio
        self.tick_track_ratio = tick_track_ratio
        self.tick_style = tick_style
        self.plot_size_thr = plot_size_thr
        self.tick_labelsize = tick_labelsize
        self._tracks: List[Track] = []

        self._check_init_values()

    def _check_init_values(self) -> None:
        """Check initial values"""
        if self.align_type not in ("left", "center", "right"):
            err_msg = f"Invalid align type '{self.align_type}'."
            raise ValueError(err_msg)

        if self.tick_style is not None and self.tick_style not in ("axis", "bar"):
            err_msg = f"Invalid tick type '{self.tick_style}'."
            raise ValueError(err_msg)

    def _get_track_offset(self, track: Track) -> int:
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
            raise NotImplementedError()

    @property
    def _track_name2offset(self) -> Dict[str, int]:
        """Track name & offset dict"""
        track_name2offset = {}
        for track in self.get_tracks(subtrack=True):
            track_name2offset[track.name] = self._get_track_offset(track)
        return track_name2offset

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
        track_names = [t.name for t in self.get_tracks()]
        feature_track_idx1 = track_names.index(feature_track_name1)
        feature_track_idx2 = track_names.index(feature_track_name2)
        if abs(feature_track_idx1 - feature_track_idx2) == 2:
            link_track_idx = int((feature_track_idx1 + feature_track_idx2) / 2)
            target_link_track = self.get_tracks()[link_track_idx]
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

    def add_feature_track(
        self,
        name: str,
        size: int,
        labelsize: int = 20,
        labelcolor: str = "black",
        labelmargin: float = 0.01,
        linewidth: int = 1,
        linecolor: str = "grey",
        link_track_ratio: Optional[float] = None,
    ) -> FeatureTrack:
        """Add feature track

        Add feature track, and also add link track between feature tracks
        if other feature tracks already exist.

        Parameters
        ----------
        name : str
            Track name
        size : int
            Track size
        labelsize : int, optional
            Track label size
        labelcolor : str, optional
            Track label color
        labelmargin : flaot, optional
            Track label margin
        linewidth : int, optional
            Track line width
        linecolor : str, optional
            Track line color
        link_track_ratio : Optional[float], optional
            Link track ratio. By default, the link_track_ratio value set
            when GenomeViz was instantiated is used.

        Returns
        -------
        feature_track : FeatureTrack
            Feature track
        """
        # Check specified track name is unique or not
        if name in [t.name for t in self.get_tracks()]:
            err_msg = f"track.name='{name}' is already exists."
            raise ValueError(err_msg)
        # Add link track between feature tracks
        if len(self.get_tracks()) > 0:
            if link_track_ratio is None:
                link_track_ratio = self.link_track_ratio
            link_track_name = f"{self._tracks[-1].name}-{name}"
            link_track = LinkTrack(link_track_name, self.track_spines, link_track_ratio)
            self._tracks.append(link_track)
        # Add feature track
        feature_track = FeatureTrack(
            name,
            size,
            labelsize,
            labelcolor,
            labelmargin,
            linewidth,
            linecolor,
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
        alpha: float = 0.8,
        v: Optional[float] = None,
        vmin: float = 0,
        vmax: float = 100,
        curve: bool = False,
        size_ratio: float = 1.0,
        patch_kws: Optional[Dict[str, Any]] = None,
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
        alpha : float, optional
            Color transparency
        v : Optional[float], optional
            Value for color interpolation
        vmin : float, optional
            Min value for color interpolation
        vmax : float, optional
            Max value for color interpolation
        curve : bool, optional
            If True, bezier curve link is plotted
        size_ratio : float, optional
            Link size ratio to track
        patch_kws : Optional[Dict[str, Any]], optional
            Optional keyword arguments to pass to link Patch object.
            See https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html
            for detailed parameters.
        """
        link_track = self._get_link_track(track_link1[0], track_link2[0])
        tracks = [t.name for t in self.get_tracks()]
        track_idx1 = tracks.index(track_link1[0])
        track_idx2 = tracks.index(track_link2[0])
        if track_idx1 < track_idx2:
            above_track_link, below_track_link = track_link1, track_link2
        else:
            above_track_link, below_track_link = track_link2, track_link1
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
                alpha,
                v,
                vmin,
                vmax,
                curve,
                size_ratio,
                patch_kws,
            )
        )

    def get_track(self, track_name: str) -> Track:
        """Get track by name

        Parameters
        ----------
        track_name : str
            Target track name

        Returns
        -------
        track : Track
            Target track
        """
        name2track = {t.name: t for t in self.get_tracks()}
        if track_name not in name2track.keys():
            err_msg = f"track.name='{track_name}' is not found."
            raise ValueError(err_msg)
        return name2track[track_name]

    def get_tracks(self, subtrack: bool = False) -> List[Track]:
        """Get tracks

        Parameters
        ----------
        subtrack : bool, optional
            If True, include feature subtracks

        Returns
        -------
        tracks : List[Track]
            Track list
        """
        tracks = []
        for track in self._tracks:
            if isinstance(track, FeatureTrack):
                tracks.append(track)
                if subtrack:
                    tracks.extend(track.subtracks)
            else:
                tracks.append(track)
        return tracks

    @property
    def top_track(self) -> FeatureTrack:
        """Top feature track"""
        feature_tracks = [t for t in self.get_tracks() if isinstance(t, FeatureTrack)]
        if len(feature_tracks) == 0:
            err_msg = "No track found. Can't access 'top_track' property."
            raise ValueError(err_msg)
        return feature_tracks[0]

    @property
    def bottom_track(self) -> FeatureTrack:
        """Bottom feature track"""
        feature_tracks = [t for t in self.get_tracks() if isinstance(t, FeatureTrack)]
        if len(feature_tracks) == 0:
            err_msg = "No track found. Can't access 'bottom_track' property."
            raise ValueError(err_msg)
        return feature_tracks[-1]

    @property
    def max_track_size(self) -> int:
        """Max track size"""
        if len(self.get_tracks()) == 0:
            err_msg = "No track found. Can't access 'max_track_size' property."
            raise ValueError(err_msg)
        return max([track.size for track in self.get_tracks()])

    def plotfig(self, dpi: int = 100) -> Figure:
        """Plot figure

        Parameters
        ----------
        dpi : int, optional
            DPI

        Returns
        -------
        figure : Figure
            Plot figure result
        """
        if len(self.get_tracks()) == 0:
            raise ValueError("No tracks are defined for plotting figure.")

        # Set tick track if required
        self._tracks = [t for t in self.get_tracks() if not isinstance(t, TickTrack)]
        max_track_size = self.max_track_size
        if self.tick_style is not None:
            self._tracks.append(
                TickTrack(
                    max_track_size,
                    self.tick_labelsize,
                    self.track_spines,
                    self.tick_track_ratio,
                    self.tick_style,
                )
            )

        # Set figure
        track_num = len(self.get_tracks(subtrack=True))
        figsize = (self.fig_width, self.fig_track_height * track_num)
        tight_layout = False if track_num < 3 else True
        figure: Figure = plt.figure(
            figsize=figsize, facecolor="white", dpi=dpi, tight_layout=tight_layout
        )

        # Set gridspec
        height_ratios = [t.ratio for t in self.get_tracks(subtrack=True)]
        spec = gridspec.GridSpec(
            nrows=track_num, ncols=1, height_ratios=height_ratios, hspace=0
        )

        # Plot each track
        plot_length_thr = max_track_size * self.plot_size_thr
        for idx, track in enumerate(self.get_tracks(subtrack=True)):
            # Create new track subplot
            xlim, ylim = (0, max_track_size), track.ylim
            ax: Axes = figure.add_subplot(
                spec[idx], xlim=xlim, ylim=ylim, fc="none", zorder=track.zorder
            )
            track_offset = self._get_track_offset(track)
            track._ax, track._offset = ax, track_offset
            # Set 'spines' and 'ticks' visibility
            for spine, display in track.spines_params.items():
                ax.spines[spine].set_visible(display)
            ax.tick_params(**track.tick_params)

            if isinstance(track, FeatureTrack):
                # Plot track label
                if track.labelsize != 0:
                    margin = -max_track_size * track.labelmargin
                    ax.text(margin, 0, **track.label_params)
                # Plot track scale line
                xmin, xmax = track_offset, track.size + track_offset
                ax.hlines(0, xmin, xmax, track.linecolor, linewidth=track.linewidth)

                for feature in [f + track_offset for f in track.features]:
                    # Don't plot too small feature (To reduce drawing time)
                    if feature.length < plot_length_thr:
                        continue
                    # Plot feature patch & label text
                    feature.plot_feature(ax, max_track_size, ylim)
                    feature.plot_label(ax, ylim)

            elif isinstance(track, FeatureSubTrack):
                # No specific plans to implementation at this time
                pass

            elif isinstance(track, LinkTrack):
                for link in track.links:
                    # Don't plot too small link (To reduce drawing time)
                    length1, length2 = link.track_length1, link.track_length2
                    if 0 < length1 < plot_length_thr or 0 < length2 < plot_length_thr:
                        continue
                    link = link.add_offset(self._track_name2offset)
                    link.plot_link(ax, ylim)

            elif isinstance(track, TickTrack):
                if self.tick_style == "axis":
                    ax.xaxis.set_major_locator(MaxNLocator(10, steps=[1, 2, 5, 10]))
                    ax.xaxis.set_major_formatter(track.tick_formatter)
                elif self.tick_style == "bar":
                    avg_height_ratio = sum(height_ratios) / len(height_ratios)
                    track.height_scale = avg_height_ratio / self.tick_track_ratio
                    ymin, ycenter, ymax = track.ymin, track.ycenter, track.ymax
                    line_kws = {"colors": "black", "linewidth": 1, "clip_on": False}
                    ax.hlines(ycenter, track.xmin, track.xmax, **line_kws)
                    ax.vlines(track.xmin, ymin, ymax, **line_kws)
                    ax.vlines(track.xmax, ymin, ymax, **line_kws)
                    ax.text(**track.scalebar_text_params)

            else:
                raise NotImplementedError()

        return figure

    def savefig(
        self,
        savefile: Union[str, Path],
        dpi: int = 100,
        pad_inches: float = 0.5,
    ) -> None:
        """Save figure to file

        Parameters
        ----------
        savefile : Union[str, Path]
            Save file
        dpi : int, optional
            DPI
        pad_inches : float, optional
            Padding inches
        """
        figure = self.plotfig(dpi=dpi)
        figure.savefig(
            fname=savefile,
            dpi=dpi,
            pad_inches=pad_inches,
            bbox_inches="tight",
        )

    def print_tracks_info(self, detail=False) -> None:
        """Print tracks info (Mainly for debugging work)

        Parameters
        ----------
        detail : bool, optional
            If True, also output feature and link details
        """
        for idx, track in enumerate(self.get_tracks(subtrack=True), 1):
            # Print track common info
            class_name = track.__class__.__name__
            print(f"\n# Track{idx:02d}: Name='{track.name}' ({class_name})")
            size, ratio, zorder = track.size, track.ratio, track.zorder
            print(f"# Size={size}, Ratio={ratio}, Zorder={zorder}", end="")

            # Print each track specific info
            if isinstance(track, FeatureTrack):
                print(f", FeatureNumber={len(track.features)}")
                if detail:
                    for feature in track.features:
                        print(feature)
            elif isinstance(track, FeatureSubTrack):
                print()
            elif isinstance(track, LinkTrack):
                print(f", LinkNumber={len(track.links)}")
                if detail:
                    for link in track.links:
                        print(link)
            elif isinstance(track, TickTrack):
                print()

    def set_colorbar(
        self,
        figure: Figure,
        bar_colors: List[str] = ["grey", "red"],
        alpha: float = 0.8,
        vmin: float = 0,
        vmax: float = 100,
        bar_height: float = 0.2,
        bar_width: float = 0.01,
        bar_left: float = 1.01,
        bar_bottom: float = 0.05,
        bar_label: str = "",
        bar_labelsize: float = 15,
        tick_labelsize: float = 10,
    ) -> None:
        """Set colorbars to figure

        Set colorbars for similarity links between genome tracks

        Parameters
        ----------
        figure : Figure
            Matplotlib figure
        bar_colors : List[str], optional
            Bar color list
        alpha : float, optional
            Color transparency
        vmin : float, optional
            Colorbar min value
        vmax : float, optional
            Colorbar max value
        bar_height : float, optional
            Colorbar height
        bar_width : float, optional
            Colorbar width
        bar_left : float, optional
            Colorbar left position
        bar_bottom : float, optional
            Colorbar bottom position
        bar_label : str, optional
            Colorbar label name
        bar_labelsize : float, optional
            Colorbar label size
        tick_labelsize : float, optional
            Colorbar tick label size
        """
        for cnt, color in enumerate(bar_colors):
            left = bar_left + bar_width * cnt
            cbar_ax = figure.add_axes([left, bar_bottom, bar_width, bar_height])

            def to_nearly_white(color: str, nearly_value: float = 0.1) -> str:
                """Convert target color to nearly white"""
                cmap = colors.LinearSegmentedColormap.from_list("m", ("white", color))
                return colors.to_hex(cmap(nearly_value))

            nearly_white = to_nearly_white(color)
            cmap = colors.LinearSegmentedColormap.from_list("m", (nearly_white, color))
            norm = colors.Normalize(vmin=vmin, vmax=vmax)
            cb_kws = {"orientation": "vertical", "ticks": []}
            cb = ColorbarBase(cbar_ax, cmap=cmap, norm=norm, alpha=alpha, **cb_kws)
            if cnt == len(bar_colors) - 1:
                ticks = [vmin, vmax]
                labels = [f"{t}%" for t in ticks]
                cb.set_ticks(ticks, labels=labels, fontsize=tick_labelsize)
                x, y = 2.0, (vmin + vmax) / 2
                text_kws = {"rotation": 90, "ha": "left", "va": "center"}
                cbar_ax.text(x, y, bar_label, size=bar_labelsize, **text_kws)
