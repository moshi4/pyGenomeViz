from __future__ import annotations

from typing import List

from pygenomeviz.feature import Feature
from pygenomeviz.link import Link


class Track:
    """Track BaseClass"""

    def __init__(
        self,
        name: str,
        size: int,
    ):
        """Track constructor

        Args:
            name (str): Track name
            size (int): Track size
        """
        self.name = name
        self.size = size


class FeatureTrack(Track):
    """FeatureTrack Class"""

    def __init__(
        self,
        name: str,
        size: int,
        labelsize: int = 30,
        linewidth: int = 2,
    ):
        """FeatureTrack constructor

        Args:
            name (str): Track name
            size (int): Track size
            labelsize (int, optional): Track label size. Defaults to 30.
            linewidth (int, optional): Track line width. Defaults to 2.
        """
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
        """Add feature to track

        Args:
            start (int): Feature start postion
            end (int): Feature end position
            strand (int): Feature strand
            label (str, optional): Feature label. Defaults to "".
            plotstyle (str, optional): Feature plot style. Defaults to "bigarrow".
            facecolor (str, optional): Feature face color. Defaults to "orange".
            edgecolor (str, optional): Feature edge color. Defaults to "black".
        """
        self.features.append(
            Feature(start, end, strand, label, plotstyle, facecolor, edgecolor)
        )


class LinkTrack(Track):
    """LinkTrack Class"""

    def __init__(
        self,
        name: str,
        size: int,
        labelsize: int = 0,
        linewidth: int = 0,
    ):
        """LinkTrack constructor

        Args:
            name (str): Track name
            size (int): Track size
            labelsize (int, optional): Track label size. Defaults to 0.
            linewidth (int, optional): Track line width. Defaults to 0.
        """
        super().__init__(name, size)
        self.labelsize = labelsize
        self.linewidth = linewidth
        self.links: List[Link] = []

    def add_link(self, link: Link) -> None:
        """Add link to track

        Args:
            link (Link): FeatureTrack to FeatureTrack Link
        """
        self.links.append(link)
