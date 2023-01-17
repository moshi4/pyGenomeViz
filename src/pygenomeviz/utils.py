from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from io import StringIO, TextIOWrapper
from pathlib import Path
from urllib.request import urlretrieve

import matplotlib.pyplot as plt
import numpy as np
from Bio import Entrez
from matplotlib.colors import to_hex

from pygenomeviz.config import DATASETS, EXAMPLE_GFF_FILES, GITHUB_DATA_URL


def load_example_gff(
    filename: str,
    cache_dir: str | Path | None = None,
    overwrite_cache: bool = False,
) -> Path:
    """Load pygenomeviz example GFF file

    Load example GFF file from https://github.com/moshi4/pygenomeviz-data/
    and cache GFF file in local directory (Default: `~/.cache/pygenomeviz/`).

    List of example GFF filename
    - `enterobacteria_phage.gff`
    - `escherichia_coli.gff.gz`

    Parameters
    ----------
    filename : str
        GFF filename (e.g. `enterobacteria_phage.gff`)
    cache_dir : str | Path | None, optional
        Output cache directory (Default: `~/.cache/pygenomeviz/`)
    overwrite_cache : bool, optional
        If True, overwrite cached GFF file

    Returns
    -------
    gff_file : Path
        GFF file
    """
    # Check specified filename exists or not
    if filename not in EXAMPLE_GFF_FILES:
        err_msg = f"{filename=} not found."
        raise ValueError(err_msg)

    # Cache local directory
    if cache_dir is None:
        package_name = __name__.split(".")[0]
        cache_base_dir = Path.home() / ".cache" / package_name
        cache_dir = cache_base_dir / "gff"
        os.makedirs(cache_dir, exist_ok=True)
    else:
        cache_dir = Path(cache_dir)

    # Download GFF file
    gff_file_url = GITHUB_DATA_URL + f"gff/{filename}"
    gff_file_path = cache_dir / filename
    if overwrite_cache or not gff_file_path.exists():
        urlretrieve(gff_file_url, gff_file_path)

    return gff_file_path


def load_dataset(
    name: str,
    cache_dir: str | Path | None = None,
    overwrite_cache: bool = False,
) -> tuple[list[Path], list[DatasetLink]]:
    """Load pygenomeviz example genbank dataset

    Load genbank datasets from https://github.com/moshi4/pygenomeviz-data
    and cache datasets in local directory (Default: `~/.cache/pygenomeviz/`).

    List of dataset name
    - `escherichia_phage`
    - `erwinia_phage`
    - `enterobacteria_phage`
    - `escherichia_coli`
    - `mycoplasma_gallisepticum`

    Parameters
    ----------
    name : str
        Dataset name (e.g. `escherichia_phage`)

    cache_dir : str | Path | None, optional
        Output cache directory (Default: `~/.cache/pygenomeviz/`)

    overwrite_cache : bool
        If True, overwrite cached dataset

    Returns
    -------
    gbk_files, links : tuple[list[Path], list[DatasetLink]]
        Genbank files, DatasetLink list
    """
    # Check specified name dataset exists or not
    if name not in DATASETS.keys():
        err_msg = f"'{name}' dataset not found."
        raise ValueError(err_msg)

    # Dataset cache local directory
    if cache_dir is None:
        package_name = __name__.split(".")[0]
        cache_base_dir = Path.home() / ".cache" / package_name
        cache_dir = cache_base_dir / name
        os.makedirs(cache_dir, exist_ok=True)
    else:
        cache_dir = Path(cache_dir)

    # Download & cache dataset
    gbk_files: list[Path] = []
    links: list[DatasetLink] = []
    for filename in DATASETS[name]:
        file_url = GITHUB_DATA_URL + f"{name}/{filename}"
        file_path = cache_dir / filename
        if overwrite_cache or not file_path.exists():
            urlretrieve(file_url, file_path)
        if file_path.suffix in (".gb", ".gbk", ".gbff"):
            gbk_files.append(file_path)
        else:
            links = DatasetLink.load(file_path)

    return gbk_files, links


@dataclass
class DatasetLink:
    """Dataset Link DataClass (Only used in `load_dataset()` function)"""

    ref_name: str
    ref_start: int
    ref_end: int
    query_name: str
    query_start: int
    query_end: int
    identity: float

    @staticmethod
    def load(link_file: Path) -> list[DatasetLink]:
        """Load pyGenomeViz dataset link file

        Parameters
        ----------
        link_file : Path
            pyGenomeViz dataset link file

        Returns
        -------
        links : list[DatasetLink]
            DatasetLink list
        """
        with open(link_file) as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader)
            links: list[DatasetLink] = []
            for row in reader:
                rname, rstart, rend = row[7], int(row[0]), int(row[1])
                qname, qstart, qend = row[8], int(row[2]), int(row[3])
                ident = float(row[6])
                links.append(
                    DatasetLink(rname, rstart, rend, qname, qstart, qend, ident)
                )
        return links


def fetch_genbank_by_accid(
    accid: str,
    gbk_outfile: str | Path | None = None,
    email: str | None = None,
) -> TextIOWrapper:
    """Fetch genbank text by 'Accession ID'

    Parameters
    ----------
    accid : str
        Accession ID
    gbk_outfile : str | Path | None, optional
        If file path is set, write fetch data to file
    email : str | None, optional
        Email address to notify download limitation (Required for bulk download)

    Returns
    -------
    TextIOWrapper
        Genbank data

    Examples
    --------
    >>> gbk_fetch_data = fetch_genbank_by_accid("JX128258.1")
    >>> gbk = Genbank(gbk_fetch_data)
    """
    Entrez.email = "" if email is None else email
    gbk_fetch_data: TextIOWrapper = Entrez.efetch(
        db="nucleotide",
        id=accid,
        rettype="gbwithparts",
        retmode="text",
    )
    if gbk_outfile is not None:
        gbk_text = gbk_fetch_data.read()
        with open(gbk_outfile, "w") as f:
            f.write(gbk_text)
        gbk_fetch_data = StringIO(gbk_text)

    return gbk_fetch_data


class ColorCycler:
    """Color Cycler Class"""

    counter = 0
    cmap = plt.get_cmap("tab10")

    def __new__(cls, n: int | None = None) -> str:
        """Get hexcolor cyclically from cmap by counter or user specified number

        `ColorCycler()` works same as `ColorCycler.get_color()` (syntactic sugar)

        Parameters
        ----------
        n : int | None, optional
            Number for color cycle. If None, counter class variable is used.

        Returns
        -------
        hexcolor : str
            Cyclic hexcolor string
        """
        return cls.get_color(n)

    @classmethod
    def reset_cycle(cls) -> None:
        """Reset cycle counter"""
        cls.counter = 0

    @classmethod
    def set_cmap(cls, name: str) -> None:
        """Set colormap (Default: `tab10`)"""
        cls.cmap = plt.get_cmap(name)
        cls.counter = 0

    @classmethod
    def get_color(cls, n: int | None = None) -> str:
        """Get hexcolor cyclically from cmap by counter or user specified number

        Parameters
        ----------
        n : int | None, optional
            Number for color cycle. If None, counter class variable is used.

        Returns
        -------
        hexcolor : str
            Cyclic hexcolor string
        """
        if n is None:
            n = cls.counter
            cls.counter += 1
        return to_hex(cls.cmap(n % cls.cmap.N), keep_alpha=True)

    @classmethod
    def get_color_list(cls, n: int | None = None) -> list[str]:
        """Get hexcolor list of colormap

        Parameters
        ----------
        n : int | None, optional
            If n is None, all(=cmap.N) hexcolors are extracted from colormap.
            If n is specified, hexcolors are extracted from n equally divided colormap.

        Returns
        -------
        hexcolor_list : list[str]
            Hexcolor list
        """
        if n is None:
            cmap_idx_list = list(range(0, cls.cmap.N))
        elif n > 0:
            cmap_idx_list = [int(i) for i in np.linspace(0, cls.cmap.N, n)]
        else:
            raise ValueError(f"{n=} is invalid number (Must be 'n > 0').")

        return [to_hex(cls.cmap(i), keep_alpha=True) for i in cmap_idx_list]
