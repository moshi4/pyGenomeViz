from __future__ import annotations

import uuid
from copy import deepcopy
from dataclasses import dataclass
from typing import Any

from matplotlib.patches import Patch

from pygenomeviz.exception import LinkRangeError
from pygenomeviz.patches import Link
from pygenomeviz.segment import FeatureSegment
from pygenomeviz.track import FeatureTrack, Track
from pygenomeviz.utils.plot import plot_patches


class LinkTrack(Track):
    """Link Track Class"""

    def __init__(
        self,
        name: str,
        *,
        ratio: float = 1.0,
        upper_feature_track: FeatureTrack,
        lower_feature_track: FeatureTrack,
    ):
        """
        Parameters
        ----------
        name : str
            Track name
        upper_feature_track : FeatureTrack
            Upper feature track
        lower_feature_track : FeatureTrack
            Lower feature track
        ratio : float, optional
            Track size ratio
        """
        super().__init__(name, ratio=ratio, zorder=0.0)

        self._upper_feature_track = upper_feature_track
        self._lower_feature_track = lower_feature_track

        self._link_record_list: list[LinkRecord] = []
        self._gid2link_dict: dict[str, dict[str, Any]] = {}

    @property
    def upper_feature_track(self) -> FeatureTrack:
        """Upper feature track"""
        return self._upper_feature_track

    @property
    def lower_feature_track(self) -> FeatureTrack:
        """Lower feature track"""
        return self._lower_feature_track

    @property
    def link_record_list(self) -> list[LinkRecord]:
        """Link record list"""
        return self._link_record_list

    @property
    def gid2link_dict(self) -> dict[str, dict[str, Any]]:
        """gid & link dict"""
        return self._gid2link_dict

    def add_link(
        self,
        upper_seg: FeatureSegment,
        upper_start: int,
        upper_end: int,
        lower_seg: FeatureSegment,
        lower_start: int,
        lower_end: int,
        *,
        v: float | None = None,
        size: float = 1.0,
        curve: bool = False,
        **kwargs,
    ) -> None:
        """Add link

        Parameters
        ----------
        upper_seg : FeatureSegment
            Upper segment
        upper_start : int
            Upper segment link start position
        upper_end : int
            Upper segment link end postion
        lower_seg : FeatureSegment
            Lower segment
        lower_start : int
            Lower segment link start position
        lower_end : int
            Lower segment link end postion
        v : float | None, optional
            Identity value for color interpolation
        size : float, optional
            Link vertical size ratio for track
        curve : bool, optional
            Curve or not
        **kwargs: dict, optional
            Patch properties (e.g. `ec="black", lw=0.5, hatch="//", ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        # Check link position is within segment range
        if not upper_seg.is_within_range((upper_start, upper_end)):
            err_msg = f"{upper_start=}, {upper_end=} is invalid ({upper_seg})"
            raise LinkRangeError(err_msg)
        if not lower_seg.is_within_range((lower_start, lower_end)):
            err_msg = f"{lower_start=}, {lower_end=} is invalid ({lower_seg})"
            raise LinkRangeError(err_msg)

        link_record = LinkRecord(
            track1=self.upper_feature_track,
            seg1=upper_seg,
            start1=upper_start,
            end1=upper_end,
            track2=self.lower_feature_track,
            seg2=lower_seg,
            start2=lower_start,
            end2=lower_end,
            v=v,
            ylim=(self.ylim[0] * size, self.ylim[1] * size),
            curve=curve,
            patch_kws=kwargs,
        )
        self._link_record_list.append(link_record)
        self._gid2link_dict[link_record.gid] = link_record.to_dict()

    def plot_links(self, fast_render: bool = True) -> None:
        """Plot links

        Parameters
        ----------
        fast_render : bool, optional
            Enable fast rendering using PatchCollection plot style.
        """
        patches: list[Patch] = [record.to_patch() for record in self.link_record_list]
        plot_patches(patches, self.ax, fast_render)


@dataclass
class LinkRecord:
    track1: FeatureTrack
    seg1: FeatureSegment
    start1: int
    end1: int
    track2: FeatureTrack
    seg2: FeatureSegment
    start2: int
    end2: int
    v: float | None = None
    ylim: tuple[float, float] = (-1, 1)
    curve: bool = False
    patch_kws: dict[str, Any] | None = None

    def __post_init__(self):
        # Generate uuid for patch group id
        self._gid = f"Link-{uuid.uuid4().hex}"

    @property
    def gid(self) -> str:
        """Group ID"""
        return self._gid

    @property
    def length1(self) -> int:
        """Length1"""
        return max(self.start1, self.end1) - min(self.start1, self.end1)

    @property
    def length2(self) -> int:
        """Length2"""
        return max(self.start2, self.end2) - min(self.start2, self.end2)

    def to_patch(self) -> Link:
        """Convert to link patch

        Returns
        -------
        link_patch : Link
            Link patch
        """
        self.patch_kws = {} if self.patch_kws is None else deepcopy(self.patch_kws)
        self.patch_kws.update(gid=self.gid)

        return Link(
            start1=self.seg1.transform_coord(self.start1),
            end1=self.seg1.transform_coord(self.end1),
            start2=self.seg2.transform_coord(self.start2),
            end2=self.seg2.transform_coord(self.end2),
            ylim=self.ylim,
            curve=self.curve,
            **self.patch_kws,
        )

    def to_dict(self) -> dict[str, Any]:
        """Convert to dict for tooltip display

        Returns
        -------
        link_dict : dict[str, Any]
            link dict
        """
        return dict(
            track1=self.track1.name,
            track2=self.track2.name,
            segment1=self.seg1.name,
            segment2=self.seg2.name,
            start1=self.start1,
            start2=self.start2,
            end1=self.end1,
            end2=self.end2,
            length1=self.length1,
            length2=self.length2,
            identity=self.v if self.v else "na",
        )
