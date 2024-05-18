from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

from pygenomeviz.typing import AlnMethod


@dataclass
class FigureConfig:
    """Figure Plot Config Class"""

    width: float = 15
    track_height: float = 1.0
    feature_track_ratio: float = 0.25
    link_track_ratio: float = 1.0
    label_size: int = 20
    range_label_size: int = 15
    track_align_type: str = "center"
    scale_style: str | None = None
    seg_space_ratio: float = 0.02


@dataclass
class FeatureConfig:
    """Feature Plot Config Class"""

    types: list[str] = field(default_factory=lambda: ["CDS"])
    type2plotstyle: dict[str, str] = field(default_factory=lambda: dict(CDS="arrow"))
    type2color: dict[str, str] = field(default_factory=lambda: dict(CDS="orange"))
    line_width: float = 0.0
    pseudo_color: str = "lightgrey"
    label_target_track: str = "top"
    label_type: str | None = None
    label_size: int = 10
    label_filter_words: list[str] = field(default_factory=list)

    @property
    def label_filter_func(self) -> Callable[[str], str]:
        """Label filter function (from label_filter_words)"""

        def label_filter(label: str) -> str:
            filter_words = self.label_filter_words + ["hypothetical"]
            for filter_word in filter_words:
                if filter_word.strip() in label:
                    return ""
            return label

        return label_filter


@dataclass
class AlignConfig:
    """Alignment & Link Config Class"""

    method: AlnMethod = None
    min_length: int = 0
    min_identity: float = 0.0
    curve: bool = False
    colorbar_height: float = 0.3
    normal_link_color: str = "grey"
    inverted_link_color: str = "red"


@dataclass
class PgvGuiPlotConfig:
    """pyGenomeViz GUI Plot Config Class"""

    fig: FigureConfig
    feat: FeatureConfig
    aln: AlignConfig
    name2seqid2range: dict[str, dict[str, tuple[int, int]]]
