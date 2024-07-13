from __future__ import annotations

import uuid
from copy import deepcopy
from typing import TYPE_CHECKING, Any, Callable, overload

import numpy as np
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation

from pygenomeviz.exception import FeatureRangeError
from pygenomeviz.typing import HPos, PlotStyle, VPos

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from pygenomeviz.track import FeatureTrack


class FeatureSegment:
    """Feature Segment Class"""

    def __init__(
        self,
        name: str,
        start: int,
        end: int,
        feature_track: FeatureTrack,
    ):
        """
        Parameters
        ----------
        name : str
            Segment name
        start : int
            Segment start position
        end : int
            Segment end position
        feature_track : FeatureTrack
            Parent feature track
        """
        self._name = name
        self._start = start
        self._end = end
        self._feature_track = feature_track

        self._features: list[SeqFeature] = []
        self._exon_features: list[SeqFeature] = []
        self._text_kws_list: list[dict[str, Any]] = []
        self._gid2feature_dict: dict[str, dict[str, Any]] = {}

    ############################################################
    # Property
    ############################################################

    @property
    def name(self) -> str:
        """Segment name"""
        return self._name

    @property
    def start(self) -> int:
        """Segment start position"""
        return self._start

    @property
    def end(self) -> int:
        """Segment end position"""
        return self._end

    @property
    def range(self) -> tuple[int, int]:
        """Segment (start, end) range"""
        return (self.start, self.end)

    @property
    def size(self) -> int:
        """Segment size"""
        return self.end - self.start

    @property
    def feature_track(self) -> FeatureTrack:
        """Parent feature track"""
        return self._feature_track

    @property
    def track_start(self) -> int:
        """Segment start position in track"""
        pos = 0
        for idx, segment in enumerate(self.feature_track.segments):
            start_pos = pos + self.feature_track.offset
            if segment.name == self.name:
                break
            if idx < len(self.feature_track.segments) - 1:
                pos += segment.size + self.feature_track.spaces[idx]
        return start_pos

    @property
    def track_end(self) -> int:
        """Segment end position in track"""
        return self.track_start + self.size

    @property
    def gid2feature_dict(self) -> dict[str, dict[str, Any]]:
        """gid & feature dict (Sort by start coordinate)"""
        return dict(sorted(self._gid2feature_dict.items(), key=lambda v: v[1]["start"]))

    @property
    def transform_features(self) -> list[SeqFeature]:
        """Coordinate transformed features

        Segment-level coordinate is transformed to track-level coordinate.
        """
        return list(map(self._transform_feature, self._features))

    @property
    def transform_exon_features(self) -> list[SeqFeature]:
        """Coordinate transformed exon features

        Segment-level coordinate is transformed to track-level coordinate.
        """
        return list(map(self._transform_feature, self._exon_features))

    @property
    def transform_text_kws_list(self) -> list[dict[str, Any]]:
        """Coordinate transformed text keywords list"""
        plot_text_kws_list: list[dict[str, Any]] = []
        for text_kws in self._text_kws_list:
            plot_text_kws = deepcopy(text_kws)
            plot_text_kws["x"] = self.transform_coord(plot_text_kws["x"])
            plot_text_kws_list.append(plot_text_kws)
        return plot_text_kws_list

    ############################################################
    # Public Method
    ############################################################

    def is_within_range(self, pos: int | tuple[int, int]) -> bool:
        """Check target pos is within segment range"""
        if isinstance(pos, int):
            return self.start <= pos <= self.end
        else:
            min_pos, max_pos = min(pos), max(pos)
            return self.start <= min_pos <= max_pos <= self.end

    @overload
    def transform_coord(self, x: int) -> int: ...
    @overload
    def transform_coord(self, x: float) -> float: ...
    @overload
    def transform_coord(self, x: NDArray) -> NDArray[np.float64]: ...

    def transform_coord(
        self, x: int | float | NDArray
    ) -> int | float | NDArray[np.float64]:
        """Transform segment-level coordinate to track-level coordinate

        Parameters
        ----------
        x : int | float| NDArray
            Segment level coordinate(s)

        Returns
        -------
        track_coord : int | float | NDArray[np.float64]
            Track level coordinate(s)
        """
        offset = self.track_start - self.start
        err_msg = f"{x=} is outside the segment range ({self})"
        if isinstance(x, (int, float)):
            if not self.start <= x <= self.end:
                raise ValueError(err_msg)
            return x + offset
        else:
            x = np.array(x)
            if np.any(x < self.start) or np.any(x > self.end):  # type: ignore
                raise ValueError(err_msg)
            return (np.array(x) + offset).astype(np.float64)

    def add_text(
        self,
        x: float,
        text: str,
        *,
        size: float = 12,
        vpos: VPos = "top",
        hpos: HPos = "left",
        ymargin: float = 0.2,
        rotation: float = 45,
        **kwargs,
    ) -> None:
        """Add text

        Parameters
        ----------
        x : float
            Text x coordinate
        text : str
            Text content
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
            Text properties (e.g. `color="red", ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.text.html>
        """
        # Ignore if text content is blank or size <= 0
        if text == "" or size <= 0:
            return

        # Check x coordinate is valid or not
        if not self.start <= x <= self.end:
            raise ValueError(f"{x=} is invalid ({self.start=}, {self.end})")

        vpos2y = dict(top=1 + ymargin, center=0, bottom=-(1 + ymargin))
        vpos2va = dict(top="bottom", center="center", bottom="top")
        hpos2ha = dict(left="left", center="center", right="right")

        text_kws = dict(
            x=x,
            y=vpos2y[vpos],
            s=text,
            size=size,
            va=vpos2va[vpos],
            ha=hpos2ha[hpos],
            rotation=rotation,
            rotation_mode="anchor",
            **kwargs,
        )
        self._text_kws_list.append(text_kws)

    def add_sublabel(
        self,
        text: str | None = None,
        *,
        size: float = 12,
        pos: str = "bottom-left",
        ymargin: float = 0.2,
        rotation: float = 0,
        **kwargs,
    ) -> None:
        """Add sublabel

        Parameters
        ----------
        text : str | None, optional
            Text content. If None, `{start:,} - {end:,} bp` label is set.
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
        # Check pos is valid or not
        vpos, hpos = pos.split("-")
        vpos_types, hpos_types = ("top", "bottom"), ("left", "center", "right")
        if vpos not in vpos_types or hpos not in hpos_types:
            err_msg = f"{pos=} is invalid pattern. "
            err_msg += "position must be '[top|bottom]-[left|center|right]'"
            raise ValueError(err_msg)

        # Set default sublabel text
        if text is None:
            text = f"{self.start:,} - {self.end:,} bp"

        hpos2x = dict(
            left=self.start,
            center=(self.start + self.end) / 2,
            right=self.end,
        )

        self.add_text(
            hpos2x[hpos],
            text,
            size=size,
            vpos=vpos,  # type: ignore
            hpos=hpos,  # type: ignore
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
        plotstyle : PlotStyle, optional
            Feature plot style (`bigarrow`|`arrow`|`bigbox`|`box`|`bigrbox`|`rbox`)
        arrow_shaft_ratio : float, optional
            Arrow shaft size ratio
        label : str, optional
            Feature label text
        text_kws : dict[str, Any] | None, optional
            `segment.add_text()` method keyword arguments
            (e.g. `dict(size=12, color="red", ...)`)
        **kwargs : dict, optional
            Patch properties (e.g. `fc="red", lw=0.5, hatch="//", ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        text_kws = {} if text_kws is None else deepcopy(text_kws)

        # Plot feature
        feature = SeqFeature(SimpleLocation(start, end, strand))
        self.add_features(
            feature,
            plotstyle=plotstyle,
            arrow_shaft_ratio=arrow_shaft_ratio,
            **kwargs,
        )

        # Plot text
        label_pos = (start + end) / 2
        self.add_text(label_pos, label, **text_kws)

    def add_features(
        self,
        features: SeqFeature | list[SeqFeature],
        *,
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
        text_kws : dict[str, Any] | None, optional
            `segment.add_text()` method keyword arguments
            (e.g. `dict(color="red", ...)`)
        **kwargs : dict, optional
            Patch properties (e.g. `fc="red", lw=0.5, hatch="//", ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        text_kws = {} if text_kws is None else deepcopy(text_kws)

        # Set default label handler
        def default_label_handler(label: str) -> str:
            return "" if "hypothetical" in label.lower() else label

        if label_handler is None:
            label_handler = default_label_handler

        if isinstance(features, SeqFeature):
            features = [features]

        for feature in features:
            try:
                # Check feature is within segment range
                self._check_feature_within_segment(feature)
            except FeatureRangeError:
                if ignore_outside_range:
                    continue
                else:
                    raise

            # Update feature qualifiers for feature patch plot
            gid = f"Feature-{uuid.uuid4().hex}"
            kwargs.update(gid=gid)
            self._add_gid2feature_dict(gid, feature, extra_tooltip)
            feature.qualifiers.update(
                plotstyle=plotstyle,
                arrow_shaft_ratio=arrow_shaft_ratio,
                patch_kws=deepcopy(kwargs),
            )
            self._features.append(feature)

            # Plot feature label
            label = feature.qualifiers.get(label_type, [""])[0]
            label = label_handler(label)
            start, end = int(feature.location.start), int(feature.location.end)  # type: ignore
            label_pos = (start + end) / 2
            self.add_text(label_pos, label, **text_kws)

    def add_exon_feature(
        self,
        locs: list[tuple[int, int]],
        strand: int = 1,
        *,
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
        text_kws = {} if text_kws is None else deepcopy(text_kws)

        # Plot exon feature
        if len(locs) == 1:
            start, end = locs[0]
            feature = SeqFeature(SimpleLocation(start, end, strand))
        else:
            feature_locs = [SimpleLocation(*loc, strand) for loc in locs]
            feature = SeqFeature(CompoundLocation(feature_locs))
        self.add_exon_features(
            feature,
            plotstyle=plotstyle,
            arrow_shaft_ratio=arrow_shaft_ratio,
            patch_kws=patch_kws,
            intron_patch_kws=intron_patch_kws,
        )

        # Plot text
        start, end = int(feature.location.start), int(feature.location.end)  # type: ignore
        label_pos = (start + end) / 2
        self.add_text(label_pos, label, **text_kws)

    def add_exon_features(
        self,
        features: SeqFeature | list[SeqFeature],
        *,
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
        """Add exon features (BioPython SeqFeature)

        Parameters
        ----------
        features : SeqFeature | list[SeqFeature]
            BioPython SeqFeature or SeqFeature list
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
        patch_kws = {} if patch_kws is None else deepcopy(patch_kws)
        intron_patch_kws = {} if intron_patch_kws is None else intron_patch_kws
        text_kws = {} if text_kws is None else deepcopy(text_kws)

        # Set default label handler
        def default_label_handler(label: str) -> str:
            return "" if "hypothetical" in label.lower() else label

        if label_handler is None:
            label_handler = default_label_handler

        if isinstance(features, SeqFeature):
            features = [features]

        for feature in features:
            try:
                # Check feature is within segment range
                self._check_feature_within_segment(feature)
            except FeatureRangeError:
                if ignore_outside_range:
                    continue
                else:
                    raise

            # Update feature qualifiers for feature patch plot
            gid = f"Feature-{uuid.uuid4().hex}"
            patch_kws.update(gid=gid)
            self._add_gid2feature_dict(gid, feature, extra_tooltip)
            feature.qualifiers.update(
                plotstyle=plotstyle,
                arrow_shaft_ratio=arrow_shaft_ratio,
                patch_kws=deepcopy(patch_kws),
                intron_patch_kws=intron_patch_kws,
            )
            self._exon_features.append(feature)

            # Plot feature label
            label = feature.qualifiers.get(label_type, [""])[0]
            label = label_handler(label)
            start, end = int(feature.location.start), int(feature.location.end)  # type: ignore
            label_pos = (start + end) / 2
            self.add_text(label_pos, label, **text_kws)

    ############################################################
    # Private Method
    ############################################################

    def _check_feature_within_segment(self, feature: SeqFeature) -> None:
        """Check if feature is within segment

        Parameters
        ----------
        feature : SeqFeature
            BioPython SeqFeature

        Raises
        ------
        FeatureOutsideRangeError
            feature is not within segment range
        """
        start, end = int(feature.location.start), int(feature.location.end)  # type: ignore
        if not self.start <= start <= end <= self.end:
            feature_location = str(feature.location)
            segment_range = f"{self.start} - {self.end}"
            err_msg = f"{feature_location=} is invalid ({segment_range=})"
            raise FeatureRangeError(err_msg)

    def _transform_feature(self, feature: SeqFeature) -> SeqFeature:
        """Transform segment-level feature coordinate to track-level

        Parameters
        ----------
        feature : SeqFeature
            BioPython SeqFeature

        Returns
        -------
        feature : SeqFeature
            Transformaed feature
        """
        locs: list[SimpleLocation] = []
        for loc in feature.location.parts:
            start = self.transform_coord(int(loc.start))  # type: ignore
            end = self.transform_coord(int(loc.end))  # type: ignore
            locs.append(SimpleLocation(start, end, loc.strand))
        transform_feature = SeqFeature(
            location=CompoundLocation(locs) if len(locs) >= 2 else locs[0],
            type=feature.type,
            qualifiers=deepcopy(feature.qualifiers),
        )
        return transform_feature

    def _add_gid2feature_dict(
        self,
        gid: str,
        feature: SeqFeature,
        extra_tooltip: dict[str, str] | None = None,
    ) -> None:
        """Add gid & feature dict

        Parameters
        ----------
        gid : str
            Group id
        feature : SeqFeature
            BioPython SeqFeature
        extra_tooltip : dict[str, str] | None, optional
            Extra tooltip dict
        """
        extra_tooltip = {} if extra_tooltip is None else deepcopy(extra_tooltip)

        start, end = int(feature.location.start), int(feature.location.end)  # type: ignore
        strand = "-" if feature.location.strand == -1 else "+"
        location = f"{start:,} - {end:,} ({strand})"

        self._gid2feature_dict[gid] = dict(
            gid=gid,
            track=self.feature_track.label,
            segment=self.name,
            start=start,
            end=end,
            strand=strand,
            location=location,
            length=end - start,
            type="na" if feature.type == "" else feature.type,
            gene=feature.qualifiers.get("gene", ["na"])[0],
            protein_id=feature.qualifiers.get("protein_id", ["na"])[0],
            product=feature.qualifiers.get("product", ["na"])[0],
            pseudo="pseudo" in feature.qualifiers or "pseudogene" in feature.qualifiers,
            translation=feature.qualifiers.get("translation", ["na"])[0],
            extra=extra_tooltip,
        )

    def __str__(self):
        seg_name = self.name
        seg_size = self.size
        seg_range = f"({self.start} - {self.end})"
        return f"{seg_name=}, {seg_size=}, {seg_range=}"

    def __repr__(self):
        return str(self)
