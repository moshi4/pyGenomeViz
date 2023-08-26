from __future__ import annotations

import time
from io import StringIO
from pathlib import Path
from typing import TYPE_CHECKING

import streamlit as st

from pygenomeviz import Genbank

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
        filename = Path(gbk_file.name).stem
        return Genbank(
            StringIO(gbk_file.getvalue().decode("utf-8")),
            name=filename.replace(" ", "_").replace("|", "_"),
        )


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
