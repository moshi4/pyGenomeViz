from __future__ import annotations

import textwrap
from collections.abc import Sequence
from copy import deepcopy
from typing import TYPE_CHECKING, Any, Callable, Mapping, overload

import numpy as np
from Bio.SeqFeature import SeqFeature
from matplotlib.patches import Patch

from pygenomeviz.exception import SegmentNotFoundError, SubTrackNotFoundError
from pygenomeviz.patches import PLOTSTYLE2PATCH, Intron
from pygenomeviz.segment import FeatureSegment
from pygenomeviz.track import Track
from pygenomeviz.typing import HPos, PlotStyle, TrackAlignType, VPos
from pygenomeviz.utils.plot import plot_patches

if TYPE_CHECKING:
    from numpy.typing import NDArray


class FeatureTrack(Track):
    """Feature Track Class"""

    def __init__(
        self,
        name: str,
        seg_name2range: Mapping[str, tuple[int, int]],
        *,
        ratio: float = 1.0,
        space: float | list[float] = 0.01,
        offset: int | TrackAlignType = "left",
        labelsize: float = 20,
        labelmargin: float = 0.01,
        align_label: bool = True,
        label_kws: dict[str, Any] | None = None,
        line_kws: dict[str, Any] | None = None,
    ):
        """
        Parameters
        ----------
        name : str
            Track name
        seg_name2range : Mapping[str, tuple[int, int]]
            Segment name & range dict
        ratio : float, optional
            Track size ratio
        space : float | list[float], optional
            Space ratio between segments
        offset : int | TrackAlignType, optional
            Offset int value or TrackAlignType (`left`|`center`|`right`)
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
        """
        label_kws = {} if label_kws is None else deepcopy(label_kws)
        line_kws = {} if line_kws is None else deepcopy(line_kws)
        line_kws.setdefault("color", "grey")
        line_kws.setdefault("lw", 1.0)

        super().__init__(name, ratio=ratio, zorder=1.0)

        segments: list[FeatureSegment] = []
        for seg_name, range in seg_name2range.items():
            segment = FeatureSegment(seg_name, *range, self)
            segments.append(segment)

        # Check space list length
        if isinstance(space, (list, tuple)):
            if len(space) != len(seg_name2range) - 1:
                raise ValueError(f"{len(space)=} is invalid!!")

        self._segments = segments
        self._space = space
        self._offset = offset
        self._labelsize = labelsize
        self._labelmargin = labelmargin
        self._align_label = align_label
        self._label_kws = label_kws
        self._line_kws = line_kws
        self._subtracks: list[FeatureSubTrack] = []

        self._label: str | None = None
        self._segment_sep_text_kws_list: Sequence[dict[str, Any] | None] = []

        self._max_track_total_seg_size: int | None = None

    ############################################################
    # Property
    ############################################################

    @property
    def label(self) -> str:
        """Track label (By default, `track.label` = `track.name`)"""
        return self.name if self._label is None else self._label

    @property
    def offset(self) -> int:
        """Track offset"""
        offset = self._offset
        if isinstance(offset, str):
            if offset == "left":
                return 0
            elif offset == "center":
                return int((max(self.xlim) - self.plot_size) / 2)
            elif offset == "right":
                return max(self.xlim) - self.plot_size
            else:
                raise ValueError(f"{offset=} is invalid!!")
        else:
            if not offset >= 0:
                raise ValueError(f"offset must be greater than 0 ({offset=}).")
            return offset

    @property
    def segments(self) -> list[FeatureSegment]:
        """Segments"""
        return self._segments

    @property
    def subtracks(self) -> list[FeatureSubTrack]:
        """Subtracks"""
        return self._subtracks

    @property
    def total_seg_size(self) -> int:
        """Total segment size"""
        return sum([seg.size for seg in self.segments])

    @property
    def spaces(self) -> list[int]:
        """Spaces between segments"""
        spaces: list[int] = []
        if isinstance(self._space, (list, tuple)):
            for space in self._space:
                if 0 <= space < 1:
                    spaces.append(int(self.max_track_total_seg_size * space))
                else:
                    spaces.append(int(space))
        else:
            for _ in range(len(self.segments) - 1):
                if 0 <= self._space < 1:
                    spaces.append(int(self.max_track_total_seg_size * self._space))
                else:
                    spaces.append(int(self._space))
        return spaces

    @property
    def max_track_total_seg_size(self) -> int:
        """Max track total segment size (Use space calculation)"""
        if self._max_track_total_seg_size is None:
            raise ValueError("'max_track_total_seg_size' is not defined!!")
        else:
            return self._max_track_total_seg_size

    @property
    def plot_size(self) -> int:
        """Plot x size (`total_seg_size` + `sum(spaces)`)"""
        return self.total_seg_size + sum(self.spaces)

    ############################################################
    # Public Method
    ############################################################

    def set_max_track_total_seg_size(self, max_track_total_seg_size: int) -> None:
        """Set max track total segment size

        This method is expected to be called within the GenomeViz instance
        to update track status. General users should not use this method.

        Parameters
        ----------
        max_track_total_seg_size : int
            Max track total segment size
        """
        self._max_track_total_seg_size = max_track_total_seg_size

    def set_label(self, label: str) -> None:
        """Set track label (By default, `track.label` = `track.name`)

        Parameters
        ----------
        label : str
            Track label
        """
        self._label = label

    def set_segment_sep(
        self,
        sep: bool | list[bool] = True,
        *,
        symbol: str = "//",
        size: float = 20,
        color: str = "grey",
        **kwargs,
    ) -> None:
        """Set segment separator symbol text

        Parameters
        ----------
        sep : bool | list[bool]
            If True, insert separator text between all segments.
            If list[bool], insert separator text between segments where True.
        symbol : str, optional
            Separator symbol text
        size : float, optional
            Separator symbol size
        color : str, optional
            Separator symbol color
        **kwargs : dict, optional
            Text properties
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.text.html>
        """
        # Check list sep length
        if isinstance(sep, (list, tuple)):
            if len(sep) != len(self.spaces):
                raise ValueError(f"{len(sep)=} is invalid!!")

        # Convert bool sep to list sep
        if isinstance(sep, bool):
            sep = [sep] * len(self.spaces)

        # Set segment separator symbol text kws
        self._segment_sep_text_kws_list = []
        text_kws = dict(s=symbol, size=size, color=color, **kwargs)
        text_kws.setdefault("va", "center")
        text_kws.setdefault("ha", "center")
        text_kws.setdefault("weight", "ultralight")
        for set_sep in sep:
            if set_sep:
                self._segment_sep_text_kws_list.append(text_kws)
            else:
                self._segment_sep_text_kws_list.append(None)

    def add_subtrack(
        self,
        name: str | None = None,
        *,
        ratio: float = 1.0,
        ylim: tuple[int, int] = (0, 100),
    ) -> FeatureSubTrack:
        """Add subtrack for user-defined plot axes

        Parameters
        ----------
        name : str
            Track name
        ratio : float, optional
            Subtrack size ratio to feature track
        ylim : tuple[int, int], optional
            Axes ylim

        Returns
        -------
        subtrack : FeatureSubTrack
            Subtrack
        """
        default_name = f"{self.name}_subtrack{len(self.subtracks) + 1:02d}"
        name = default_name if name is None else name

        # Check track name duplication
        subtrack_names = [t.name for t in self.subtracks]
        if name in subtrack_names:
            raise ValueError(f"{name=} subtrack is already exists!!")

        subtrack = FeatureSubTrack(name, ratio=self.ratio * ratio, feature_track=self)
        subtrack.set_xlim(self.xlim)
        subtrack.set_ylim(ylim)
        self._subtracks.append(subtrack)

        return subtrack

    def get_subtrack(self, name: str | None = None) -> FeatureSubTrack:
        """Get subtrack by name

        If no subtrack found, raise error.

        Parameters
        ----------
        name : str | None, optional
            Target subtrack name. If None, first subtrack is returned.

        Returns
        -------
        subtrack : FeatureSubTrack
            Target subtrack
        """
        if len(self.subtracks) == 0:
            raise SubTrackNotFoundError("Failed to get subtrack. No subtrack found.")
        if name is None:
            return self.subtracks[0]
        else:
            name2subtrack = {t.name: t for t in self.subtracks}
            if name not in name2subtrack:
                raise SubTrackNotFoundError(f"{name=} subtrack not found.")
            return name2subtrack[name]

    def get_segment(
        self,
        name: str | None = None,
    ) -> FeatureSegment:
        """Get segment by name

        Parameters
        ----------
        name : str | None
            Target segment name. If None, first segment is returned.

        Returns
        -------
        segment : FeatureSegment
            Target segment
        """
        if name is None:
            return self.segments[0]
        else:
            name2segment = {seg.name: seg for seg in self.segments}
            if name not in name2segment:
                err_msg = f"{name=} segment not found (track_name='{self.name}')."
                raise SegmentNotFoundError(err_msg)
            return name2segment[name]

    def add_text(
        self,
        x: float,
        text: str,
        *,
        target_seg: str | None = None,
        size: float = 15,
        vpos: VPos = "top",
        hpos: HPos = "left",
        ymargin: float = 0.2,
        rotation: float = 45,
        **kwargs,
    ) -> None:
        """Add text to track segment

        Parameters
        ----------
        x : float
            Text x coordinate
        text : str
            Text content
        target_seg : str | None, optional
            Target segment name. If None, first segment is selected.
        size : float, optional
            Text size
        vpos : str, optional
            Vertical position (`top`|`center`|`bottom`)
        hpos : str, optional
            Horizontal position (`left`|`center`|`right`)
        ymargin : float, optional
            Y margin
        rotation : float, optional
            Text rotation
        **kwargs : dict, optional
            `segment.add_text()` method keyword arguments (e.g. `color="red", ...`)
        """
        segment = self.get_segment(target_seg)
        segment.add_text(
            x,
            text,
            size=size,
            vpos=vpos,
            hpos=hpos,
            ymargin=ymargin,
            rotation=rotation,
            **kwargs,
        )

    def add_sublabel(
        self,
        text: str | None = None,
        *,
        target_seg: str | None = None,
        size: float = 12,
        pos: str = "bottom-left",
        ymargin: float = 0.2,
        rotation: float = 0,
        **kwargs,
    ) -> None:
        """Add sublabel to corners of the track segment

        Parameters
        ----------
        text : str | None, optional
            Text content
        target_seg : str | None, optional
            Target segment name. If None, first segment is selected.
        size : float, optional
            Text size
        pos : str, optional
            Label position ([`top`|`bottom`]-[`left`|`center`|`right`])
        ymargin : float, optional
            Y margin
        rotation : float, optional
            Text rotation
        **kwargs : dict, optional
            `segment.add_text()` method keyword arguments (e.g. `color="red", ...`)
        """
        segment = self.get_segment(target_seg)
        segment.add_sublabel(
            text,
            size=size,
            pos=pos,
            ymargin=ymargin,
            rotation=rotation,
            **kwargs,
        )

    def add_feature(
        self,
        start: int,
        end: int,
        strand: int = 1,
        *,
        target_seg: str | None = None,
        plotstyle: PlotStyle = "arrow",
        arrow_shaft_ratio: float = 0.5,
        label: str = "",
        text_kws: dict[str, Any] | None = None,
        **kwargs,
    ) -> None:
        """Add feature

        Parameters
        ----------
        start : int
            Start position
        end : int
            End position
        strand : int, optional
            Feature strand
        target_seg : str | None, optional
            Target segment name. If None, first segment is selected.
        plotstyle : PlotStyle, optional
            Feature plot style (`bigarrow`|`arrow`|`bigbox`|`box`|`bigrbox`|`rbox`)
        arrow_shaft_ratio : float, optional
            Arrow shaft size ratio
        label : str, optional
            Feature label
        text_kws : dict[str, Any] | None, optional
            `segment.add_text()` method keyword arguments
            (e.g. `dict(size=12, color="red", ...)`)
        **kwargs : dict, optional
            Patch properties (e.g. `fc="red", lw=0.5, hatch="//", ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        segment = self.get_segment(target_seg)
        segment.add_feature(
            start,
            end,
            strand,
            plotstyle=plotstyle,
            arrow_shaft_ratio=arrow_shaft_ratio,
            label=label,
            text_kws=text_kws,
            **kwargs,
        )

    def add_features(
        self,
        features: SeqFeature | list[SeqFeature],
        *,
        target_seg: str | None = None,
        plotstyle: PlotStyle = "arrow",
        arrow_shaft_ratio: float = 0.5,
        label_type: str | None = None,
        label_handler: Callable[[str], str] | None = None,
        extra_tooltip: dict[str, str] | None = None,
        ignore_outside_range: bool = False,
        text_kws: dict[str, Any] | None = None,
        **kwargs,
    ) -> None:
        """Add features (BioPython SeqFeature)

        Parameters
        ----------
        features : SeqFeature | list[SeqFeature]
            BioPython SeqFeature or SeqFeature list
        target_seg : str | None, optional
            Target segment name. If None, first segment is selected.
        plotstyle : PlotStyle, optional
            Feature plot style (`bigarrow`|`arrow`|`bigbox`|`box`|`bigrbox`|`rbox`)
        arrow_shaft_ratio : float, optional
            Arrow shaft size ratio
        label_type : str | None, optional
            Label type (e.g. `gene`,`protein_id`,`product`,etc...)
        label_handler : Callable[[str], str] | None, optional
            Label handler function to customize label display.
            If None, set label handler to exclude labels containing `hypothetical`.
        extra_tooltip : dict[str, str] | None, optional
            Extra tooltip dict for html figure
        ignore_outside_range : bool, optional
            If True and the feature position is outside the range of the track segment,
            ignore it without raising an error.
        text_kws : dict[str, Any] | None, optional
            `segment.add_text()` method keyword arguments
            (e.g. `dict(size=12, color="red", ...)`)
        **kwargs : dict, optional
            Patch properties (e.g. `fc="red", lw=0.5, hatch="//", ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        segment = self.get_segment(target_seg)
        segment.add_features(
            features,
            plotstyle=plotstyle,
            arrow_shaft_ratio=arrow_shaft_ratio,
            label_type=label_type,
            label_handler=label_handler,
            extra_tooltip=extra_tooltip,
            ignore_outside_range=ignore_outside_range,
            text_kws=text_kws,
            **kwargs,
        )

    def add_exon_feature(
        self,
        locs: list[tuple[int, int]],
        strand: int = 1,
        *,
        target_seg: str | None = None,
        plotstyle: PlotStyle = "arrow",
        arrow_shaft_ratio: float = 0.5,
        label: str = "",
        patch_kws: dict[str, Any] | None = None,
        intron_patch_kws: dict[str, Any] | None = None,
        text_kws: dict[str, Any] | None = None,
    ) -> None:
        """Add exon feature

        Parameters
        ----------
        locs : list[tuple[int, int]]
            Exon locations (e.g. `[(0, 100), (200, 300), (350, 400)]`)
        strand : int, optional
            Feature strand
        target_seg : str | None, optional
            Target segment name. If None, first segment is selected.
        plotstyle : PlotStyle, optional
            Feature plot style (`bigarrow`|`arrow`|`bigbox`|`box`|`bigrbox`|`rbox`)
        arrow_shaft_ratio : float, optional
            Arrow shaft size ratio
        label : str, optional
            Feature label
        patch_kws : dict[str, Any] | None, optional
            Exon patch properties (e.g. `dict(fc="red", lw=0.5, hatch="//", ...)`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        intron_patch_kws : dict[str, Any] | None, optional
            Intron patch properties (e.g. `dict(color="red", lw=2.0, ...)`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        text_kws : dict[str, Any] | None, optional
            `segment.add_text()` method keyword arguments
            (e.g. `dict(size=12, color="red", ...)`)
        """
        segment = self.get_segment(target_seg)
        segment.add_exon_feature(
            locs,
            strand,
            plotstyle=plotstyle,
            arrow_shaft_ratio=arrow_shaft_ratio,
            label=label,
            patch_kws=patch_kws,
            intron_patch_kws=intron_patch_kws,
            text_kws=text_kws,
        )

    def add_exon_features(
        self,
        features: SeqFeature | list[SeqFeature],
        *,
        target_seg: str | None = None,
        plotstyle: PlotStyle = "arrow",
        arrow_shaft_ratio: float = 0.5,
        label_type: str | None = None,
        label_handler: Callable[[str], str] | None = None,
        extra_tooltip: dict[str, str] | None = None,
        ignore_outside_range: bool = False,
        patch_kws: dict[str, Any] | None = None,
        intron_patch_kws: dict[str, Any] | None = None,
        text_kws: dict[str, Any] | None = None,
    ) -> None:
        """Add exon features

        Parameters
        ----------
        features : SeqFeature | list[SeqFeature]
            BioPython SeqFeature or SeqFeature list
        target_seg : str | None, optional
            Target segment name. If None, first segment is selected.
        plotstyle : PlotStyle, optional
            Feature plot style (`bigarrow`|`arrow`|`bigbox`|`box`|`bigrbox`|`rbox`)
        arrow_shaft_ratio : float, optional
            Arrow shaft size ratio
        label_type : str | None, optional
            Label type (e.g. `gene`,`protein_id`,`product`, etc...)
        label_handler : Callable[[str], str] | None, optional
            Label handler function to customize label display.
            If None, set label handler to exclude labels containing `hypothetical`.
        extra_tooltip : dict[str, str] | None, optional
            Extra tooltip dict for html figure
        ignore_outside_range : bool, optional
            If True and the feature position is outside the range of the track segment,
            ignore it without raising an error.
        patch_kws : dict[str, Any] | None, optional
            Exon patch properties (e.g. `dict(fc="red", lw=0.5, hatch="//", ...)`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        intron_patch_kws : dict[str, Any] | None, optional
            Intron patch properties (e.g. `dict(color="red", lw=2.0, ...)`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        text_kws : dict[str, Any] | None, optional
            `segment.add_text()` method keyword arguments
            (e.g. `dict(size=12, color="red", ...)`)
        """
        segment = self.get_segment(target_seg)
        segment.add_exon_features(
            features,
            plotstyle=plotstyle,
            arrow_shaft_ratio=arrow_shaft_ratio,
            label_type=label_type,
            label_handler=label_handler,
            extra_tooltip=extra_tooltip,
            ignore_outside_range=ignore_outside_range,
            patch_kws=patch_kws,
            intron_patch_kws=intron_patch_kws,
            text_kws=text_kws,
        )

    @overload
    def transform_coord(self, x: int, *, target_seg: str | None = None) -> int: ...
    @overload
    def transform_coord(self, x: float, *, target_seg: str | None = None) -> float: ...
    @overload
    def transform_coord(
        self, x: NDArray, *, target_seg: str | None = None
    ) -> NDArray[np.float64]: ...

    def transform_coord(
        self,
        x: int | float | NDArray,
        *,
        target_seg: str | None = None,
    ) -> int | float | NDArray[np.float64]:
        """Transform segment-level coordinate to track-level coordinate

        Parameters
        ----------
        x : int | float | NDArray
            Segment-level coordinate(s)
        target_seg : str | None, optional
            Target segment name. If None, first segment is selected.

        Returns
        -------
        transform_x : int | float| NDArray[np.float64]
            Track-level coordinate(s)
        """
        seg = self.get_segment(target_seg)
        return seg.transform_coord(x)

    def plot_all(self, fast_render: bool = True) -> None:
        """Plot all objects (Expected to be called in `gv.plotfig()`)

        1. Plot track label
        2. Plot segment lines
        3. Plot segment separator
        4. Plot features
        5. Plot texts

        Parameters
        ----------
        fast_render : bool, optional
            Enable fast rendering using PatchCollection plot style.
        """
        self._plot_track_label()
        self._plot_segment_lines()
        self._plot_segment_sep()
        self._plot_features(fast_render)
        self._plot_exon_features(fast_render)
        self._plot_texts()

    ############################################################
    # Private Method
    ############################################################

    def _plot_track_label(self) -> None:
        """Plot track label"""
        # Calculate track label position
        if self._align_label:
            x, y = -self._labelmargin, 0.5
        else:
            first_seg_start_x = self.segments[0].track_start / self.xlim[1]
            x, y = first_seg_start_x - self._labelmargin, 0.5

        # Plot track label
        self._label_kws.update(
            ha="right",
            va="center",
            fontsize=self._labelsize,
            transform=self.ax.transAxes,
        )
        self.ax.text(x, y, self.label, **self._label_kws)

    def _plot_segment_lines(self) -> None:
        """Plot lines for each segment"""
        for seg in self.segments:
            x, y = (seg.track_start, seg.track_end), (0, 0)
            self.ax.plot(x, y, **self._line_kws)

    def _plot_segment_sep(self) -> None:
        """Plot break symbol for each segment"""
        if len(self._segment_sep_text_kws_list) == 0:
            return
        pos = 0
        for idx, space in enumerate(self.spaces):
            target_pos = pos + self.segments[idx].size + (space / 2)
            seg_sep_text_kws = self._segment_sep_text_kws_list[idx]
            if seg_sep_text_kws is not None:
                self.ax.text(target_pos, 0, **seg_sep_text_kws)
            pos += self.segments[idx].size + space

    def _plot_features(
        self,
        fast_render: bool = True,
    ) -> None:
        """Plot features for each segment

        Parameters
        ----------
        fast_render : bool, optional
            Enable fast rendering using PatchCollection plot style.
        """
        # Collect feature patches
        patches: list[Patch] = []
        for seg in self.segments:
            for f in seg.transform_features:
                start = int(f.location.parts[0].start)  # type: ignore
                end = int(f.location.parts[-1].end)  # type: ignore
                strand = int(f.location.strand)  # type: ignore
                plotstyle = str(f.qualifiers["plotstyle"])
                arrow_shaft_ratio = float(f.qualifiers["arrow_shaft_ratio"])
                patch_kws = dict(f.qualifiers["patch_kws"])
                if "arrow" in plotstyle or "rbox" in plotstyle:
                    patch_kws.update(max_size=self.xlim[1])
                if "arrow" in plotstyle:
                    patch_kws.update(shaft_ratio=arrow_shaft_ratio)

                PlotPatch = PLOTSTYLE2PATCH[plotstyle]
                patches.append(PlotPatch(start, end, strand, **patch_kws))

        plot_patches(patches, self.ax, fast_render)

    def _plot_exon_features(
        self,
        fast_render: bool = True,
    ) -> None:
        """Plot exon features for each segment

        Parameters
        ----------
        fast_render : bool, optional
            Enable fast rendering using PatchCollection plot style.
        """
        # Collect feature patches
        patches: list[Patch] = []
        for seg in self.segments:
            for f in seg.transform_exon_features:
                exon_locs, intron_locs = self._extract_exon_intron_locs(f)
                plotstyle = str(f.qualifiers["plotstyle"])
                arrow_shaft_ratio = float(f.qualifiers["arrow_shaft_ratio"])
                patch_kws = dict(f.qualifiers["patch_kws"])

                # Plot exon patches
                strand = int(f.location.strand)  # type: ignore
                exon_locs = exon_locs[::-1] if strand == -1 else exon_locs
                for idx, exon_loc in enumerate(exon_locs, 1):
                    exon_start, exon_end = exon_loc
                    if "arrow" in plotstyle or "rbox" in plotstyle:
                        patch_kws.update(max_size=self.xlim[1])
                    if "arrow" in plotstyle:
                        patch_kws.update(shaft_ratio=arrow_shaft_ratio)
                        if idx == len(exon_locs):
                            patch_kws.update(show_head=True)
                        else:
                            patch_kws.update(show_head=False)
                    Patch = PLOTSTYLE2PATCH[plotstyle]
                    patches.append(Patch(exon_start, exon_end, strand, **patch_kws))

                # Plot intron patches
                intron_patch_kws = dict(f.qualifiers["intron_patch_kws"])
                for intron_loc in intron_locs:
                    intron_start, intron_end = intron_loc
                    bigstyle = "big" in plotstyle
                    patches.append(
                        Intron(
                            intron_start,
                            intron_end,
                            strand,
                            bigstyle=bigstyle,
                            **intron_patch_kws,
                        )
                    )

        plot_patches(patches, self.ax, fast_render)

    def _plot_texts(self) -> None:
        """Plot texts"""
        for seg in self.segments:
            for text_kws in seg.transform_text_kws_list:
                self.ax.text(**text_kws)

    def _extract_exon_intron_locs(
        self,
        feature: SeqFeature,
    ) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
        """Extract exon & intron locations

        Parameters
        ----------
        feature : SeqFeature
            Exon Feature

        Returns
        -------
        exon_locs : list[tuple[int, int]]
            Exon locations
        intron_locs : list[tuple[int, int]]
            Intron locations
        """
        exon_locs: list[tuple[int, int]] = []
        intron_locs: list[tuple[int, int]] = []
        # Extract exon locations
        for loc in feature.location.parts:
            exon_start, exon_end = int(loc.start), int(loc.end)  # type: ignore
            exon_locs.append((exon_start, exon_end))
        # Extract intron locations
        for i in range(len(exon_locs) - 1):
            intron_start, intron_end = exon_locs[i][1], exon_locs[i + 1][0]
            intron_locs.append((intron_start, intron_end))

        return exon_locs, intron_locs

    def __str__(self):
        track_segments = {seg.name: (seg.start, seg.end) for seg in self.segments}
        return textwrap.dedent(
            f"""
            track_name='{self.name}' ({len(self.segments)} segments)
            {track_segments=}
            """
        )[1:-1]

    def __repr__(self):
        return str(self)


class FeatureSubTrack(Track):
    """Feature SubTrack Class"""

    def __init__(
        self,
        name: str,
        *,
        ratio: float,
        feature_track: FeatureTrack,
    ):
        """
        Parameters
        ----------
        name : str
            Track name
        ratio : float, optional
            Track size ratio
        feature_track : FeatureTrack
            Parent feature track to which subtrack belongs
        """
        super().__init__(name, ratio=ratio)
        self._feature_track = feature_track
        self.transform_coord = self.feature_track.transform_coord

    @property
    def feature_track(self) -> FeatureTrack:
        """Parent feature track to which subtrack belongs"""
        return self._feature_track

    def set_ylim(self, ylim: tuple[float, float]) -> None:
        """Set track ylim"""
        self._ylim = ylim
