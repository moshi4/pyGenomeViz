from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable


@dataclass
class FigureConfig:
    """Figure Plot Config Class"""

    width: float = 15
    track_height: float = 1.0
    tick_style: str | None = None
    align_type: str = "center"
    track_ratio: float = 0.3
    label_size: int = 20
    range_label_size: int = 15


@dataclass
class FeatureConfig:
    """Feature Plot Config Class"""

    types: list[str]
    type2color: dict[str, str]
    type2plotstyle: dict[str, str]
    label_type: str | None = None
    label_size: int = 10
    label_filter_words: list[str] = field(default_factory=list)
    show_only_top_label: bool = False

    @property
    def label_filter_func(self) -> Callable[[str], str]:
        """Label filter function (from label_filter_words)"""

        def label_filter(label: str) -> str:
            for filter_word in self.label_filter_words:
                if filter_word.strip() in label:
                    return ""
            return label

        return label_filter


@dataclass
class AlignConfig:
    """Alignment & Link Config Class"""

    method: str | None
    normal_link_color: str
    inverted_link_color: str
    curve: bool = False
    colorbar_height: float = 0.3
    min_length: int = 0
    min_identity: float = 0.0


@dataclass
class PgvConfig:
    """pyGenomeViz Plot Config Class"""

    fig: FigureConfig
    feat: FeatureConfig
    aln: AlignConfig
