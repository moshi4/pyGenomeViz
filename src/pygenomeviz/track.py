from __future__ import annotations

from typing import List, Optional

from pygenomeviz.feature import Feature
from pygenomeviz.link import Link


class Track:
    def __init__(
        self,
        name: str,
        size: int,
    ):
        """Constructor"""
        self.name = name
        self.size = size


class FeatureTrack(Track):
    """Feature Track"""

    def __init__(
        self,
        name: str,
        size: int,
        labelsize: int = 30,
        linewidth: int = 2,
    ):
        """Constructor"""
        super().__init__(name, size)
        self.labelsize = labelsize
        self.linewidth = linewidth
        self.features: List[Feature] = []

    def add_feature(
        self,
        start: int,
        end: int,
        strand: int,
        label: str = "",
        plotstyle: str = "bigarrow",
        facecolor: str = "orange",
        edgecolor: str = "black",
    ) -> None:
        """Add feature"""
        self.features.append(
            Feature(start, end, strand, label, plotstyle, facecolor, edgecolor)
        )


class LinkTrack(Track):
    """Link Track"""

    def __init__(
        self,
        name: str,
        size: int,
        labelsize: int = 0,
        linewidth: int = 0,
    ):
        """Constructor"""
        super().__init__(name, size)
        self.labelsize = labelsize
        self.linewidth = linewidth
        self.links: List[Link] = []

    def add_link(self, link: Link) -> None:
        """Add link"""
        self.links.append(link)
