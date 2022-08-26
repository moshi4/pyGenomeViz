from dataclasses import dataclass

from typing_extensions import Literal

DATASETS = {
    "escherichia_phage": [
        "JX128258.gbk",
        "MH051335.gbk",
        "MK373776.gbk",
        "MH816966.gbk",
        "link.tsv",
    ],
    "erwinia_phage": [
        "MT939486.gbk",
        "MT939487.gbk",
        "MT939488.gbk",
        "LT960552.gbk",
        "link.tsv",
    ],
    "enterobacteria_phage": [
        "NC_019724.gbk",
        "NC_024783.gbk",
        "NC_016566.gbk",
        "NC_013600.gbk",
        "NC_031081.gbk",
        "NC_028901.gbk",
        "link.tsv",
    ],
    "escherichia_coli": [
        "NC_000913.gbk",
        "NC_002695.gbk",
        "NC_011751.gbk",
        "NC_011750.gbk",
        "link.tsv",
    ],
    "mycoplasma_gallisepticum": [
        "NC_004829.gbk",
        "NC_018407.gbk",
        "NC_018408.gbk",
        "NC_018409.gbk",
        "NC_017502.gbk",
        "NC_017503.gbk",
        "link.tsv",
    ],
}


@dataclass
class LiteralTypes:
    PLOTSTYLE = Literal["bigarrow", "arrow", "bigbox", "box", "bigrbox", "rbox"]
    STRAND = Literal[1, -1, 0]
    LABELVPOS = Literal["top", "center", "bottom", "strand"]
    LABELHPOS = Literal["left", "center", "right"]
    LABELHA = Literal["left", "center", "right"]
    ALIGN_TYPE = Literal["left", "center", "right"]
    TICK_STYLE = Literal["axis", "bar", None]
