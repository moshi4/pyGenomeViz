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
        labelsize: int = 0,
        linewidth: int = 0,
    ):
        """Track constructor

        Args:
            name (str): Track name
            size (int): Track size
            labelsize (int, optional): Track label size
            linewidth (int, optional): Track line width
        """
        self.name = name
        self.size = size
        self.labelsize = labelsize
        self.linewidth = linewidth


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
            labelsize (int, optional): Track label size
            linewidth (int, optional): Track line width
        """
        super().__init__(name, size, labelsize, linewidth)
        self.features: List[Feature] = []

    def add_feature(
        self,
        start: int,
        end: int,
        strand: int,
        label: str = "",
        labelsize: int = 15,
        labelcolor: str = "black",
        plotstyle: str = "bigarrow",
        facecolor: str = "orange",
        edgecolor: str = "black",
    ) -> None:
        """Add feature to track

        Args:
            start (int): Feature start postion
            end (int): Feature end position
            strand (int): Feature strand
            label (str, optional): Feature label
            labelsize (int, optional): Feature label size
            labelcolor (int, optional): Feature label color
            plotstyle (str, optional): Feature plot style
            facecolor (str, optional): Feature face color
            edgecolor (str, optional): Feature edge color
        """
        self.features.append(
            Feature(
                start,
                end,
                strand,
                label,
                labelsize,
                labelcolor,
                plotstyle,
                facecolor,
                edgecolor,
            )
        )


class LinkTrack(Track):
    """LinkTrack Class"""

    def __init__(
        self,
        name: str,
        size: int = 0,
        labelsize: int = 0,
        linewidth: int = 0,
    ):
        """LinkTrack constructor

        Args:
            name (str): Track name
            size (int, optional): Track size
            labelsize (int, optional): Track label size
            linewidth (int, optional): Track line width
        """
        super().__init__(name, size, labelsize, linewidth)
        self.links: List[Link] = []

    def add_link(self, link: Link) -> None:
        """Add link to track

        Args:
            link (Link): FeatureTrack to FeatureTrack Link
        """
        self.links.append(link)
