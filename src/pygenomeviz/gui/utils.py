from __future__ import annotations

import gzip
import io
import time
from collections import defaultdict
from io import StringIO
from pathlib import Path
from typing import TYPE_CHECKING

import streamlit as st
from Bio.SeqFeature import SeqFeature
from matplotlib.figure import Figure

from pygenomeviz import GenomeViz
from pygenomeviz.parser import Genbank

if TYPE_CHECKING:
    from streamlit.runtime.uploaded_file_manager import UploadedFile


@st.cache_data(ttl=3600)
def load_gbk_file(gbk_file: str | Path | UploadedFile) -> Genbank:
    """Load genbank file

    Parameters
    ----------
    gbk_file : str | Path | UploadedFile
        Genbank file

    Returns
    -------
    gbk : Genbank
        Genbank parse object
    """
    if isinstance(gbk_file, (str, Path)):
        return Genbank(gbk_file)
    else:
        gbk_file_path = Path(gbk_file.name)
        trans_table = str.maketrans({" ": "_", "|": "_", "(": "[", ")": "]"})
        gbk_name = gbk_file_path.stem.translate(trans_table)

        if gbk_file_path.suffix == ".gz":
            for remove_suffix in (".gbff", ".gbk", ".gb"):
                gbk_name = gbk_name.replace(remove_suffix, "")
            with gzip.open(gbk_file, "rt", encoding="utf-8") as f:
                return Genbank(StringIO(f.read()), name=gbk_name)
        else:
            return Genbank(StringIO(gbk_file.getvalue().decode("utf-8")), name=gbk_name)


def is_st_cloud() -> bool:
    """Is launch on streamlit cloud"""
    try:
        return st.context.url.startswith("https://pygenomeviz.streamlit.app")
    except Exception:
        return False


def remove_old_files(target_dir: Path, ttl: int = 3600) -> None:
    """Remove old file in target directory

    Parameters
    ----------
    target_dir : Path
        Target directory path
    ttl : int, optional
        Time to live
    """
    for file in target_dir.glob("*"):
        elapsed_time = time.time() - file.stat().st_mtime
        if elapsed_time > ttl:
            file.unlink(missing_ok=True)


def extract_all_feature_types(gbk_list: list[Genbank], sort: bool = True) -> list[str]:
    """Extract all feature types from genbank list

    Parameters
    ----------
    gbk_list : list[Genbank]
        Genbank object list
    sort : bool, optional
        Sort feature types (`CDS`, `rRNA`, `tRNA`, ...)

    Returns
    -------
    all_feature_types : list[str]
        All feature types
    """
    # Extract all feature types set
    all_feature_types_set: set[str] = set()
    for gbk in gbk_list:
        for features in gbk.get_seqid2features(feature_type=None).values():
            for feature in features:
                all_feature_types_set.add(feature.type)
    all_feature_types = list(all_feature_types_set)
    # Sort by `CDS`, `rRNA`, `tRNA`, ...
    if sort:
        top_feature_types = []
        for feature_type in ["CDS", "rRNA", "tRNA"]:
            if feature_type in all_feature_types:
                top_feature_types.append(feature_type)
        all_feature_types = list(dict.fromkeys(top_feature_types + all_feature_types))
    return all_feature_types


def get_features_count_label(features: list[SeqFeature]) -> str:
    """Get features count label

    Parameters
    ----------
    features : list[SeqFeature]
        Target features

    Returns
    -------
    label : str
        Count label
    """
    # Count each feature type
    feature_type2count: dict[str, int] = defaultdict(int)
    for feature in features:
        feature_type2count[feature.type] += 1
    # Sort by count
    sorted_feature_type2count = {}
    for feature_type, count in sorted(
        feature_type2count.items(), key=lambda v: v[1], reverse=True
    ):
        sorted_feature_type2count[feature_type] = count
    # Create features count label
    label = "**Feature Count**  \n"
    for feature_type, count in sorted_feature_type2count.items():
        label += f"{feature_type}: {count}  \n"
    return label


def get_fig_bytes(gv: GenomeViz, fig: Figure, format: str) -> io.BytesIO:
    """Get figure bytes"""
    fig_bytes = io.BytesIO()
    if format in ("png", "svg"):
        fig.savefig(fig_bytes, format=format)
    elif format == "html":
        gv.savefig_html(fig_bytes)
    else:
        raise ValueError(f"{format=} is invalid.")
    return fig_bytes
