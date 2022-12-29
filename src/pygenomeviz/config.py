from dataclasses import dataclass
from pathlib import Path
from typing import Literal

###########################################################
# GitHub Genbank Dataset & GFF Files Config
###########################################################

GITHUB_DATA_URL = "https://raw.githubusercontent.com/moshi4/pygenomeviz-data/master/"

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

EXAMPLE_GFF_FILES = [
    "enterobacteria_phage.gff",
    "escherichia_coli.gff.gz",
]

###########################################################
# Literal Types Config
###########################################################


@dataclass
class LiteralTypes:
    PLOTSTYLE = Literal["bigarrow", "arrow", "bigbox", "box", "bigrbox", "rbox"]
    STRAND = int
    LABELVPOS = Literal["top", "center", "bottom", "strand"]
    LABELHPOS = Literal["left", "center", "right"]
    LABELHA = Literal["left", "center", "right"]
    ALIGN_TYPE = Literal["left", "center", "right"]
    TICK_STYLE = Literal["axis", "bar", None]


###########################################################
# HTML Viewer Template & Assets Config
###########################################################

_viewer_dir = Path(__file__).parent / "viewer"
_assets_dir = _viewer_dir / "assets"
_assets_files = [
    "lib/spectrum.min.css",
    "lib/jquery-ui.min.css",
    "lib/jquery.min.js",
    "lib/spectrum.min.js",
    "lib/jquery-ui.min.js",
    "lib/panzoom.min.js",
    "pgv-viewer.js",
]
ASSETS_FILES = [_assets_dir / f for f in _assets_files]
TEMPLATE_HTML_FILE = _viewer_dir / "pgv-viewer-template.html"
