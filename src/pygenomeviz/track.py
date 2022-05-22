from __future__ import annotations

from typing import List

from pygenomeviz.feature import Feature


class Track:
    def __init__(
        self,
        name: str,
        size: int,
        labelsize: int = 30,
        linewidth: int = 2,
        type: str = "feature",
    ):
        """Constructor"""
        self.name = name
        self.size = size
        self.labelsize = labelsize
        self.linewidth = linewidth
        self.type = type
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
