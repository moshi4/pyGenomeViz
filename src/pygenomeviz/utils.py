from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from io import TextIOWrapper
from pathlib import Path
from typing import List, Optional, Tuple, Union
from urllib.request import urlretrieve

import matplotlib.pyplot as plt
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
}


def load_dataset(
    name: str,
    cache_dir: Optional[Union[str, Path]] = None,
    overwrite_cache: bool = False,
) -> Tuple[List[Path], List[DatasetLink]]:
    """Load pygenomeviz example dataset

    Download and load datasets from https://github.com/moshi4/pygenomeviz-data
    and cache datasets in local directory ('~/.cache/pygenomeviz/').

    List of dataset name
    - `escherichia_phage`
    - `erwinia_phage`
    - `enterobacteria_phage`
    - `escherichia_coli`

    Parameters
    ----------
    name : str
        Dataset name (e.g. `escherichia_phage`)

    cache_dir : Optional[Union[str, Path]], optional
        Cache directory (Default: `~/.cache/pygenomeviz/`)

    overwrite_cache : bool
        If True, overwrite cached dataset

    Returns
    -------
    gbk_files, links : Tuple[List[Path], List[DatasetLink]]
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
    gbk_files: List[Path] = []
    links: List[DatasetLink] = []
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
    def load(link_file: Path) -> List[DatasetLink]:
        """Load pyGenomeViz dataset link file

        Parameters
        ----------
        link_file : Path
            pyGenomeViz dataset link file

        Returns
        -------
        links : List[DatasetLink]
            DatasetLink list
        """
        with open(link_file) as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader)
            links: List[DatasetLink] = []
            for row in reader:
                rname, rstart, rend = row[7], int(row[0]), int(row[1])
                qname, qstart, qend = row[8], int(row[2]), int(row[3])
                ident = float(row[6])
                links.append(
                    DatasetLink(rname, rstart, rend, qname, qstart, qend, ident)
                )
        return links


def fetch_genbank_by_accid(accid: str, email: Optional[str] = None) -> TextIOWrapper:
    """Fetch genbank text by 'Accession ID'

    Parameters
    ----------
    accid : str
        Accession ID
    email : str, optional
        Email address to notify download limitation (Required for bulk download)

    Returns
    -------
    TextIOWrapper
        Genbank data

    Examples
    --------
    >>> gbk_fetch_data = download_genbank_from_accid("JX128258.1")
    >>> gbk = Genbank(gbk_fetch_data)
    """
    Entrez.email = "" if email is None else email
    return Entrez.efetch(
        db="nucleotide",
        id=accid,
        rettype="gbwithparts",
        retmode="text",
    )


class ColorCycler:
    """Color Cycler Class"""

    counter = 0
    cmap = plt.get_cmap("tab10")

    def __new__(cls, n: Optional[int] = None) -> str:
        """Get hexcolor cyclically from cmap by counter or user specified number

        Parameters
        ----------
        n : Optional[int], optional
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
