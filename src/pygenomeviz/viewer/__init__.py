from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import pygenomeviz


def _concat_target_files_contents(files: list[Path], target_ext: str) -> str:
    """Concatenate target extension files contents"""
    contents = "\n"
    target_files = [file for file in files if file.suffix == target_ext]
    for target_file in target_files:
        with open(target_file) as f:
            contents += f.read() + "\n"
    return contents


_viewer_dir = Path(__file__).parent
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
_assets_files = [_assets_dir / f for f in _assets_files]

TEMPLATE_HTML_FILE = _viewer_dir / "pgv-viewer-template.html"
CSS_CONTENTS = _concat_target_files_contents(_assets_files, ".css")
JS_CONTENTS = _concat_target_files_contents(_assets_files, ".js")


def setup_viewer_html(
    svg_figure: str,
    gid2feature_tooltip: dict[str, str],
    gid2link_tooltip: dict[str, str],
) -> str:
    """Setup viewer html (Embed SVG figure, CSS & JS assets)

    Parameters
    ----------
    svg_figure : str
        SVG figure strings
    gid2feature_tooltip : dict[str, str]
        GID(Group ID) & feature tooltip dict
    gid2link_tooltip : dict[str, str]
        GID(Group ID) & link tooltip dict

    Returns
    -------
    viewer_html : str
        Viewer html strings
    """
    # Read template html file
    with open(TEMPLATE_HTML_FILE) as f:
        viewer_html = f.read()
    # Replace template strings
    return (
        viewer_html.replace("$PGV_SVG_FIG", f"\n{svg_figure}")
        .replace("$VERSION", pygenomeviz.__version__)
        .replace("$DATETIME_NOW", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        .replace("$CSS_CONTENTS", CSS_CONTENTS)
        .replace(
            "$JS_CONTENTS",
            JS_CONTENTS.replace(
                "FEATURE_TOOLTIP_JSON = {}",
                f"FEATURE_TOOLTIP_JSON = {json.dumps(gid2feature_tooltip, indent=4)}",
            ).replace(
                "LINK_TOOLTIP_JSON = {}",
                f"LINK_TOOLTIP_JSON = {json.dumps(gid2link_tooltip, indent=4)}",
            ),
        )
    )
