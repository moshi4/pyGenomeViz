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


def load_dataset(
    name: str,
    cache_dir: str | Path | None = None,
    overwrite_cache: bool = False,
) -> tuple[list[Path], list[DatasetLink]]:
    """Load pygenomeviz example dataset

    Download and load datasets from https://github.com/moshi4/pygenomeviz-data
    and cache datasets in local directory (Default: '~/.cache/pygenomeviz/').

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
        Cache directory (Default: `~/.cache/pygenomeviz/`)

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

    # Dataset GitHub URL
    base_url = "https://raw.githubusercontent.com/moshi4/pygenomeviz-data/master/"
    target_url = base_url + f"{name}/"

    # Download & cache dataset
    gbk_files: list[Path] = []
    links: list[DatasetLink] = []
    for filename in DATASETS[name]:
        file_url = target_url + filename
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
    """Dataset Link DataClass (Only used in load_dataset() function)"""

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
    def reset_cycle(cls) -> None:
        """Reset cycle counter"""
        cls.counter = 0

    @classmethod
    def set_cmap(cls, name: str) -> None:
        """Set colormap (Default: `tab10`)"""
        cls.cmap = plt.get_cmap(name)
        cls.counter = 0

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
            raise ValueError(f"n={n} is invalid number (Must be 'n > 0').")

        return [to_hex(cls.cmap(i), keep_alpha=True) for i in cmap_idx_list]
