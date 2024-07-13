from __future__ import annotations

import io
from collections.abc import Mapping, Sequence
from copy import deepcopy
from pathlib import Path
from typing import Any, Callable

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colorbar import Colorbar
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_hex
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

from pygenomeviz.exception import (
    FeatureTrackNotFoundError,
    LinkRangeError,
    LinkTrackNotFoundError,
)
from pygenomeviz.track import FeatureSubTrack, FeatureTrack, LinkTrack, Track
from pygenomeviz.typing import TrackAlignType, Unit
from pygenomeviz.utils.helper import interpolate_color, size_label_formatter
from pygenomeviz.viewer import setup_viewer_html


class GenomeViz:
    """Genome Visualization Class"""

    # By default, after saving a figure using the `savefig()` method, figure object is
    # automatically deleted to avoid memory leaks (no display on jupyter notebook)
    # If you want to display the figure on jupyter notebook using `savefig()` method,
    # set clear_savefig=False.
    clear_savefig: bool = True

    def __init__(
        self,
        *,
        fig_width: float = 15,
        fig_track_height: float = 1.0,
        track_align_type: TrackAlignType = "left",
        feature_track_ratio: float = 0.25,
        link_track_ratio: float = 1.0,
        show_axis: bool = False,
    ):
        """
        Parameters
        ----------
        fig_width : float, optional
            Figure width
        fig_track_height : float, optional
            Figure height = `fig_track_height * track number`
        track_align_type : TrackAlignType, optional
            Figure track alignment type (`left`|`center`|`right`)
        feature_track_ratio : float, optional
            Feature track size ratio
        link_track_ratio : float, optional
            Link track size ratio
        show_axis : bool, optional
            Show axis for debug purpose
        """
        self._fig_width = fig_width
        self._fig_track_height = fig_track_height
        self._track_align_type: TrackAlignType = track_align_type
        self._feature_track_ratio = feature_track_ratio
        self._link_track_ratio = link_track_ratio
        self._show_axis = show_axis

        self._tracks: list[Track] = []
        self._plot_colorbar: Callable[[Figure], None] | None = None
        self._plot_scale_bar: Callable[[Axes], None] | None = None
        self._plot_axis_ticks: Callable[[Axes], None] | None = None

    ############################################################
    # Property
    ############################################################

    @property
    def figsize(self) -> tuple[float, float]:
        """Figure size (Width, Height)"""
        track_num = len(self.get_tracks(subtrack=True))
        return (self._fig_width, self._fig_track_height * track_num)

    @property
    def feature_tracks(self) -> list[FeatureTrack]:
        """Feature tracks"""
        return [t for t in self.get_tracks() if isinstance(t, FeatureTrack)]

    @property
    def link_tracks(self) -> list[LinkTrack]:
        """Link tracks"""
        return [t for t in self.get_tracks() if isinstance(t, LinkTrack)]

    ############################################################
    # Public Method
    ############################################################

    def get_tracks(self, *, subtrack: bool = True) -> list[Track]:
        """Get tracks

        Parameters
        ----------
        subtrack : bool, optional
            If True, include subtracks in FeatureTrack

        Returns
        -------
        tracks : list[Track]
            Tracks [`FeatureTrack`|`FeatureSubTrack`|`LinkTrack`]
        """
        tracks = []
        for track in self._tracks:
            tracks.append(track)
            if subtrack and isinstance(track, FeatureTrack):
                tracks.extend(track.subtracks)
        return tracks

    def add_feature_track(
        self,
        name: str,
        segments: int
        | tuple[int, int]
        | Sequence[int | tuple[int, int]]
        | Mapping[str, int | tuple[int, int]],
        *,
        space: float | list[float] = 0.02,
        offset: int | TrackAlignType | None = None,
        labelsize: float = 20,
        labelmargin: float = 0.01,
        align_label: bool = True,
        label_kws: dict[str, Any] | None = None,
        line_kws: dict[str, Any] | None = None,
    ) -> FeatureTrack:
        """Add feature track

        Add feature track, and also add link track between feature tracks
        if other feature tracks already exist.

        Parameters
        ----------
        name : str
            Track name
        segments : int | tuple[int, int] | Sequence[int | tuple[int, int]] | Mapping[str, int | tuple[int, int]]
            Track segments definition. Segment sizes or ranges can be specified.
        space : float | list[float], optional
            Space ratio between segments.
            If `float`, all spaces are set to the same value.
            If `list[float]`, each space is set to the corresponding value
            (list size must be `len(segments) - 1`)
        offset : int | TrackAlignType | None, optional
            Offset int value or TrackAlignType (`left`|`center`|`right`)
            If None, offset is defined by GenomeViz `track_align_type` argument at initialization.
        labelsize : float, optional
            Track label size
        labelmargin : float, optional
            Track label margin
        align_label : bool, optional
            If True, align track label to the most left position.
            If False, set track label to first segment start position.
        label_kws : dict[str, Any] | None, optional
            Text properties (e.g. `dict(size=25, color="red", ...)`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.text.html>
        line_kws : dict[str, Any] | None, optional
            Axes.plot properties (e.g. `dict(color="grey", lw=0.5, ls="--", ...)`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.plot.html>

        Returns
        -------
        feature_track : FeatureTrack
            Feature track
        """  # noqa: E501
        label_kws = {} if label_kws is None else deepcopy(label_kws)
        line_kws = {} if line_kws is None else deepcopy(line_kws)

        # Check track name duplication
        if name in [t.name for t in self.get_tracks()]:
            raise ValueError(f"{name=} track is already exists!!")

        # Instantiate feature track
        feature_track = FeatureTrack(
            name,
            self._to_seg_name2range(segments),
            ratio=self._feature_track_ratio,
            space=space,
            offset=self._track_align_type if offset is None else offset,
            labelsize=labelsize,
            labelmargin=labelmargin,
            align_label=align_label,
            label_kws=label_kws,
            line_kws=line_kws,
        )

        # Add link track between feature tracks
        feature_track_num = len(self.feature_tracks)
        if feature_track_num > 0:
            upper_feature_track = self.feature_tracks[feature_track_num - 1]
            link_track_name = f"{upper_feature_track.name}-{name}"
            link_track = LinkTrack(
                link_track_name,
                ratio=self._link_track_ratio,
                upper_feature_track=upper_feature_track,
                lower_feature_track=feature_track,
            )
            self._tracks.append(link_track)

        self._tracks.append(feature_track)
        self._update_track_status()

        return feature_track

    def add_link(
        self,
        target1: tuple[str, int, int] | tuple[str, str | None, int, int],
        target2: tuple[str, int, int] | tuple[str, str | None, int, int],
        color: str = "grey",
        inverted_color: str | None = None,
        alpha: float = 0.8,
        v: float | None = None,
        vmin: float = 0,
        vmax: float = 100,
        size: float = 1.0,
        curve: bool = False,
        filter_length: int = 0,
        ignore_outside_range: bool = False,
        v_tooltip: float | None = None,
        **kwargs,
    ) -> None:
        """Add link patch to link track between adjacent feature tracks

        Parameters
        ----------
        target1 : tuple[str, int, int] | tuple[str, int | str | None, int, int]
            Target link1 `(track_name, start, end)` or `(track_name, target_segment, start, end)`
        target2 : tuple[str, int, int] | tuple[str, int | str | None, int, int]
            Target link2 `(track_name, start, end)` or `(track_name, target_segment, start, end)`
        color : str, optional
            Link color
        inverted_color : str | None, optional
            Inverted link color. If None, `color` is set.
        alpha : float, optional
            Color transparency
        v : float | None, optional
            Identity value for color interpolation. If None, no color interpolation is done.
        vmin : float, optional
            Min value for color interpolation
        vmax : float, optional
            Max value for color interpolation
        size : float, optional
            Link vertical size ratio for track
        curve : bool, optional
            If True, bezier curve link is plotted
        filter_length : int, optional
            If link length is shorter than `filter_length`, ignore it.
        ignore_outside_range : bool, optional
            If True and the link position is outside the range of the target track,
            ignore it without raising an error.
        v_tooltip: float | None, optional
            Identity value for only tooltip display.
            If no color interpolation is required, use this option instead of `v`.
        **kwargs: dict, optional
            Patch properties (e.g. `ec="black", lw=0.5, hatch="//", ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """  # noqa: E501
        # Get target link track
        if len(target1) == 3:
            target1 = (target1[0], None, target1[1], target1[2])
        if len(target2) == 3:
            target2 = (target2[0], None, target2[1], target2[2])
        link_track = self._get_target_link_track(target1, target2)

        # Get upper & lower target feature track, segment info
        track_name2target = {target1[0]: target1, target2[0]: target2}
        upper_target = track_name2target[link_track.upper_feature_track.name]
        lower_target = track_name2target[link_track.lower_feature_track.name]
        upper_seg_name, upper_start, upper_end = upper_target[1:]
        lower_seg_name, lower_start, lower_end = lower_target[1:]
        upper_seg = link_track.upper_feature_track.get_segment(upper_seg_name)
        lower_seg = link_track.lower_feature_track.get_segment(lower_seg_name)

        # Filter by link length
        upper_size = max(upper_start, upper_end) - min(upper_start, upper_end)
        lower_size = max(lower_start, lower_end) - min(lower_start, lower_end)
        if upper_size < filter_length or lower_size < filter_length:
            return

        # If target is inverted, set inverted_color to color
        is_inverted = (upper_end - upper_start) * (lower_end - lower_start) < 0
        if inverted_color is not None and is_inverted:
            color = inverted_color

        # Interpolate color if v exists
        if v is not None:
            color = interpolate_color(color, v, vmin, vmax)

        # Set patch properties
        kwargs.update(color=color, alpha=alpha)

        try:
            link_track.add_link(
                upper_seg,
                upper_start,
                upper_end,
                lower_seg,
                lower_start,
                lower_end,
                v=v if v is not None else v_tooltip,
                size=size,
                curve=curve,
                **kwargs,
            )
        except LinkRangeError:
            if ignore_outside_range:
                return
            else:
                raise

    def set_scale_bar(
        self,
        *,
        ymargin: float = 1.0,
        labelsize: float = 15,
        scale_size_label: tuple[int, str] | None = None,
    ) -> None:
        """Set scale bar

        Parameters
        ----------
        ymargin : float, optional
            Scale bar y margin
        labelsize : float, optional
            Label size
        scale_size_label : tuple[int, str] | None, optional
            Scale bar size & label tuple (e.g. `(1000, "1.0 kb")`)
            If None, scale bar size & label are automatically set.
        """

        def plot_scale_bar(lowest_track_ax: Axes):
            if scale_size_label is None:
                scale_size = lowest_track_ax.get_xticks()[1]
                scale_label = size_label_formatter(scale_size)
            else:
                scale_size, scale_label = scale_size_label
            scale_bar = AnchoredSizeBar(
                lowest_track_ax.transData,
                size=scale_size,
                label=scale_label,
                loc="upper right",
                sep=5,
                frameon=False,
                fontproperties=FontProperties(size=labelsize),
                bbox_to_anchor=(1, -ymargin),
                bbox_transform=lowest_track_ax.transAxes,
            )
            lowest_track_ax.add_artist(scale_bar)

        self._plot_scale_bar = plot_scale_bar

    def set_scale_xticks(
        self,
        *,
        ymargin: float = 1.0,
        labelsize: float = 15,
        start: int = 0,
        unit: Unit | None = None,
    ) -> None:
        """Set scale xticks

        Parameters
        ----------
        ymargin : float, optional
            X ticks y margin
        labelsize : float, optional
            Label size
        start : int, optional
            X ticks start position
        unit : Unit | None, optional
            Display unit (`Gb`|`Mb`|`Kb`|`bp`)
        """

        def plot_axis_ticks(lowest_track_ax: Axes) -> None:
            # Create new axes for axis xticks
            ticks_ax = lowest_track_ax.twiny()

            # Setup ticks properties of new axes
            x_size = max(lowest_track_ax.get_xlim()) - min(lowest_track_ax.get_xlim())
            xlim = (start, start + x_size)
            ticks_ax.set_xlim(*xlim)
            for pos in ("top", "left", "right"):
                ticks_ax.spines[pos].set_visible(False)
            ticks_ax.tick_params(top=False, bottom=False, left=False, right=False)
            ticks_ax.tick_params(labelsize=labelsize)
            ticks_ax.xaxis.set_ticks_position("bottom")
            ticks_ax.spines["bottom"].set_position(("axes", -ymargin))

            # Plot axis xticks
            xticks: list[float] = ticks_ax.get_xticks()  # type: ignore
            xticks = list(filter(lambda x: min(xlim) <= x <= max(xlim), xticks))
            xticklabels = size_label_formatter(xticks, unit)
            ticks_ax.set_xticks(xticks)
            ticks_ax.set_xticklabels(xticklabels)

            # Remove lowest track unnecessary ticks
            lowest_track_ax.set_xticks([])

        self._plot_axis_ticks = plot_axis_ticks

    def set_colorbar(
        self,
        colors: list[str] | None = None,
        *,
        alpha: float = 0.8,
        vmin: float = 0,
        vmax: float = 100,
        bar_height: float = 0.2,
        bar_width: float = 0.01,
        bar_left: float = 1.02,
        bar_bottom: float = 0,
        bar_label: str = "",
        bar_labelsize: float = 15,
        tick_labelsize: float = 10,
    ) -> None:
        """Set colorbar

        Parameters
        ----------
        colors : list[str] | None, optional
            Colors for bar
        alpha : float, optional
            Color transparency
        vmin : float, optional
            Colorbar min value
        vmax : float, optional
            Colorbar max value
        bar_height : float, optional
            Colorbar height ratio
        bar_width : float, optional
            Colorbar width ratio
        bar_left : float, optional
            Colorbar left position
        bar_bottom : float, optional
            Colorbar bottom position
        bar_label : str, optional
            Colorbar label
        bar_labelsize : float, optional
            Colorbar label size
        tick_labelsize : float, optional
            Colorbar tick label size
        """
        if bar_height <= 0 or bar_width <= 0:
            return

        # Set colors & remove duplicated colors
        colors = ["grey"] if colors is None else colors
        colors = list(dict.fromkeys([to_hex(c) for c in colors]))

        def plot_colorbar(fig: Figure) -> None:
            for cnt, color in enumerate(colors):
                # Add new axes for colorbar
                left = bar_left + bar_width * cnt
                cbar_ax = fig.add_axes((left, bar_bottom, bar_width, bar_height))
                # Set colorbar
                nealy_white = interpolate_color(color, v=0)
                cmap = LinearSegmentedColormap.from_list("m", [nealy_white, color])
                norm = Normalize(vmin=vmin, vmax=vmax)
                cb_kws = dict(orientation="vertical", ticks=[])
                cb = Colorbar(cbar_ax, cmap=cmap, norm=norm, alpha=alpha, **cb_kws)  # type: ignore
                if cnt == len(colors) - 1:
                    # Set vmin, vmax ticks label
                    ticks = [vmin, vmax]
                    labels = [f"{t}%" for t in ticks]
                    cb.set_ticks(ticks, labels=labels, fontsize=tick_labelsize)
                    # Set colorbar label
                    x, y = 2.0, (vmin + vmax) / 2
                    text_kws = dict(rotation=90, ha="left", va="center")
                    cbar_ax.text(x, y, bar_label, size=bar_labelsize, **text_kws)  # type: ignore

        self._plot_colorbar = plot_colorbar

    def plotfig(
        self,
        *,
        dpi: int = 100,
        fast_render: bool = True,
    ) -> Figure:
        """Plot figure

        Parameters
        ----------
        dpi : int, optional
            DPI
        fast_render : bool, optional
            Enable fast rendering mode using PatchCollection.
            Set fast_render=True by default, and set it to False
            when used in the `savefig_html()` method.
            Fast rendering mode cannot generate tooltips for html display.

        Returns
        -------
        fig : Figure
            Plot figure result
        """
        # Check track num
        tracks = self.get_tracks(subtrack=True)
        if len(tracks) == 0:
            raise ValueError("Failed to plot figure. No track found!!")

        # Setup figure & gridspece
        fig = plt.figure(figsize=self.figsize, dpi=dpi, facecolor="white")
        fig.tight_layout()
        height_ratios = [t.ratio for t in tracks]
        gs = GridSpec(nrows=len(tracks), ncols=1, height_ratios=height_ratios)
        gs.update(left=0, right=1, bottom=0, top=1, hspace=0, wspace=0)

        for idx, track in enumerate(tracks):
            # Create axes & set axes to track
            ax: Axes = fig.add_subplot(gs[idx])
            track.set_ax(ax, self._show_axis)

            if isinstance(track, FeatureTrack):
                track.plot_all(fast_render)
            elif isinstance(track, FeatureSubTrack):
                pass
            elif isinstance(track, LinkTrack):
                track.plot_links(fast_render)
            else:
                track_class = track.__class__.__name__
                raise NotImplementedError(f"{track_class=} is invalid track class!!")

        lowest_track_ax = tracks[-1].ax
        if self._plot_scale_bar:
            self._plot_scale_bar(lowest_track_ax)

        if self._plot_axis_ticks:
            self._plot_axis_ticks(lowest_track_ax)

        if self._plot_colorbar:
            self._plot_colorbar(fig)

        return fig

    def savefig(
        self,
        savefile: str | Path,
        *,
        dpi: int = 100,
        pad_inches: float = 0.5,
    ) -> None:
        """Save figure to file

        Parameters
        ----------
        savefile : str | Path
            Save file
        dpi : int, optional
            DPI
        pad_inches : float, optional
            Padding inches

        Warnings
        --------
        To plot a figure that settings a user-defined legend, subtracks, or annotations,
        call `fig.savefig()` instead of `gv.savefig()`.
        """
        fig = self.plotfig(dpi=dpi)
        fig.savefig(
            fname=str(savefile),
            dpi=dpi,
            pad_inches=pad_inches,
            bbox_inches="tight",
        )
        # Clear & close figure to suppress memory leak
        if self.clear_savefig:
            fig.clear()
            plt.close(fig)

    def savefig_html(
        self,
        html_outfile: str | Path | io.StringIO | io.BytesIO,
        figure: Figure | None = None,
    ) -> None:
        """Save figure in html format

        Parameters
        ----------
        html_outfile : str | Path | StringIO | BytesIO
            Output HTML file (*.html)
        figure : Figure | None, optional
            Save HTML viewer file using user customized figure.
            Set to output figure including user-specified legend, subtracks, etc.
            Target figure must be generated by `gv.plotfig(fast_render=False)`.
        """
        # Load SVG contents
        fig = self.plotfig(fast_render=False) if figure is None else figure
        svg_fig_bytes = io.BytesIO()
        fig.savefig(fname=svg_fig_bytes, format="svg")
        svg_fig_bytes.seek(0)
        svg_fig_contents = svg_fig_bytes.read().decode("utf-8")

        # Clear & close figure to suppress memory leak
        if figure is None:
            fig.clear()
            plt.close(fig)

        # Setup viewer html SVG & embed CSS, JS assets
        viewer_html = setup_viewer_html(
            svg_fig_contents,
            self._get_gid2feature_dict(),
            self._get_gid2link_dict(),
        )

        # Write viewer html contents
        if isinstance(html_outfile, io.StringIO):
            html_outfile.write(viewer_html)
        elif isinstance(html_outfile, io.BytesIO):
            html_outfile.write(bytes(viewer_html, encoding="utf-8"))
        else:
            with open(html_outfile, "w") as f:
                f.write(viewer_html)

    ############################################################
    # Private Method
    ############################################################

    def _update_track_status(self) -> None:
        """Update track status

        This method is called at the end of `add_feature_track()`
        """
        # Set `max_track_total_seg_size` & `xlim` for each feature track
        max_track_total_seg_size = max([t.total_seg_size for t in self.feature_tracks])
        for t in self.feature_tracks:
            t.set_max_track_total_seg_size(max_track_total_seg_size)

        plot_size_list = []
        for t in self.feature_tracks:
            offset = t._offset if isinstance(t._offset, int) else 0
            plot_size_list.append(t.plot_size + offset)
        for t in self.get_tracks(subtrack=True):
            t.set_xlim((0, max(plot_size_list)))

    def _to_seg_name2range(
        self,
        segments: int
        | tuple[int, int]
        | Sequence[int | tuple[int, int]]
        | Mapping[str, int | tuple[int, int]],
    ) -> dict[str, tuple[int, int]]:
        """Convert segments type to `segment name` & `range` dict"""
        # Convert segments to dict type
        if isinstance(segments, int):
            segments = [segments]
        if isinstance(segments, tuple) and len(segments) == 2:
            if isinstance(segments[0], int) and isinstance(segments[1], int):
                segments = [segments]  # type: ignore
        if isinstance(segments, (list, tuple)):
            segments = {f"seg{idx}": seg for idx, seg in enumerate(segments, 1)}  # type: ignore

        seg_name2range: dict[str, tuple[int, int]] = {}
        for seg_name, value in segments.items():  # type: ignore
            start, end = (0, value) if isinstance(value, int) else value
            if not start < end:
                raise ValueError(f"end must be larger than start ({start=}, {end=})")
            seg_name2range[seg_name] = (start, end)

        return seg_name2range

    def _get_target_link_track(
        self,
        target1: tuple[str, str | None, int, int],
        target2: tuple[str, str | None, int, int],
    ) -> LinkTrack:
        for target in (target1, target2):
            track_name, target_seg, _, _ = target
            # Check track name
            track_name2feature_track = {t.name: t for t in self.feature_tracks}
            if track_name not in track_name2feature_track:
                raise FeatureTrackNotFoundError(f"{track_name=} is not found!!")
            # Check segment name
            target_feature_track = track_name2feature_track[track_name]
            target_feature_track.get_segment(target_seg)

        target_link_track = None
        name2target = {t[0]: t for t in (target1, target2)}
        for link_track in self.link_tracks:
            upper_track_name = link_track.upper_feature_track.name
            lower_track_name = link_track.lower_feature_track.name
            if upper_track_name in name2target and lower_track_name in name2target:
                target_link_track = link_track

        if target_link_track is None:
            err_msg = f"Failed to get link track.\n{target1=}\n{target2=})\n"
            err_msg += "Target feature tracks must be adjacent feature tracks!!"
            raise LinkTrackNotFoundError(err_msg)

        return target_link_track

    def _get_gid2feature_dict(self) -> dict[str, dict[str, Any]]:
        """Get group ID & feature dict

        Returns
        -------
        gid2feature_dict : dict[str, dict[str, Any]]
            Group ID & feature dict
        """
        gid2feature_dict = {}
        for feature_track in self.feature_tracks:
            for seg in feature_track.segments:
                gid2feature_dict.update(seg.gid2feature_dict)
        return gid2feature_dict

    def _get_gid2link_dict(self) -> dict[str, dict[str, Any]]:
        """Get group ID & link dict

        Returns
        -------
        gid2link_dict : dict[str, dict[str, Any]]
            Group ID & link dict
        """
        gid2link_dict = {}
        for link_track in self.link_tracks:
            gid2link_dict.update(link_track.gid2link_dict)
        return gid2link_dict

    def __str__(self):
        ret_val = ""
        for feature_track in self.feature_tracks:
            ret_val += f"{feature_track}\n"
            for idx, seg in enumerate(feature_track.segments, 1):
                ret_val += f"  {idx}: {seg}\n"
        return ret_val
