from dataclasses import dataclass

from typing_extensions import Literal


@dataclass
class LiteralTypes:
    PLOTSTYLE = Literal["bigarrow", "arrow", "bigbox", "box", "bigrbox", "rbox"]
    STRAND = Literal[1, -1, 0]
    LABELVPOS = Literal["top", "center", "bottom", "strand"]
    LABELHPOS = Literal["left", "center", "right"]
    LABELHA = Literal["left", "center", "right"]
    ALIGN_TYPE = Literal["left", "center", "right"]
    TICK_STYLE = Literal["axis", "bar", None]
