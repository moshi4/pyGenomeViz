from __future__ import annotations

import logging
import os
from io import StringIO, TextIOWrapper
from pathlib import Path
from urllib.request import urlretrieve

from Bio import Entrez

from pygenomeviz.parser import Genbank
from pygenomeviz.typing import GenbankDatasetName, GffExampleFileName

GITHUB_DATA_URL = "https://raw.githubusercontent.com/moshi4/pygenomeviz-data-v1/main/"

GBK_DATASET = {
    "acinetobacter_phage": [
        "NC_049491.gbk",
        "NC_049492.gbk",
        "NC_049493.gbk",
        "NC_049494.gbk",
    ],
    "yersinia_phage": [
        "NC_070914.gbk",
        "NC_070915.gbk",
        "NC_070916.gbk",
        "NC_070918.gbk",
    ],
    "enterobacteria_phage": [
        "NC_013600.gbk",
        "NC_016566.gbk",
        "NC_019724.gbk",
        "NC_024783.gbk",
        "NC_028901.gbk",
        "NC_031081.gbk",
    ],
    "mycoplasma_mycoides": [
        "GCF_000023685.1.gbff",
        "GCF_000800785.1.gbff",
        "GCF_000959055.1.gbff",
        "GCF_000959065.1.gbff",
    ],
    "escherichia_coli": [
        "NC_000913.gbk.gz",
        "NC_002695.gbk.gz",
        "NC_011751.gbk.gz",
        "NC_011750.gbk.gz",
    ],
    "saccharomyces": [
        "Saccharomyces_cerevisiae.gbff.gz",
        "Saccharomyces_kudriavzevii.gbff.gz",
        "Saccharomyces_mikatae.gbff.gz",
    ],
}

GFF_FILES = [
    "enterobacteria_phage.gff",
    "mycoplasma_mycoides.gff",
    "escherichia_coli.gff.gz",
    "saccharomyces_cerevisiae.gff.gz",
]


def load_example_fasta_dataset(
    name: GenbankDatasetName,
    *,
    cache_dir: str | Path | None = None,
    overwrite_cache: bool = False,
) -> list[Path]:
    """Load pygenomeviz example fasta dataset

    Load genbank datasets from <https://github.com/moshi4/pygenomeviz-data-v1>
    and convert genbank to fasta format.
    Cache datasets in local directory (Default: `~/.cache/pygenomeviz/`).

    List of dataset name

    - `acinetobacter_phage` (4 species)
    - `yersinia_phage` (4 species)
    - `enterobacteria_phage` (6 species)
    - `mycoplasma_mycoides` (4 species)
    - `escherichia_coli` (4 species, gzip compressed)
    - `saccharomyces` (3 species, gzip compressed)

    Parameters
    ----------
    name : str
        Dataset name (e.g. `enterobacteria_phage`)
    cache_dir : str | Path | None, optional
        Output cache directory (Default: `~/.cache/pygenomeviz/`)
    overwrite_cache : bool, optional
        If True, overwrite cached dataset
    quiet : bool, optional
        If True, no print log on screen.

    Returns
    -------
    fasta_files : list[Path]
        Fasta files
    """
    logger = logging.getLogger(__name__)

    # Download genbank dataset
    gbk_files = load_example_genbank_dataset(name)
    gbk_list = list(map(Genbank, gbk_files))

    # Dataset cache local directory
    if cache_dir is None:
        package_name = __name__.split(".")[0]
        cache_base_dir = Path.home() / ".cache" / package_name
        cache_dir = cache_base_dir / "fasta" / name
        os.makedirs(cache_dir, exist_ok=True)
    else:
        cache_dir = Path(cache_dir)

    # Convert genbank to fasta format
    fasta_files: list[Path] = []
    for gbk in gbk_list:
        fasta_file = cache_dir / f"{gbk.name}.fna"
        if overwrite_cache or not fasta_file.exists():
            gbk.write_genome_fasta(fasta_file)
            logger.info(f"Save fasta file in '{fasta_file}'")
        else:
            logger.info(f"Cached fasta files found in '{fasta_file}'")
        fasta_files.append(fasta_file)

    return fasta_files


def load_example_genbank_dataset(
    name: GenbankDatasetName,
    *,
    cache_dir: str | Path | None = None,
    overwrite_cache: bool = False,
) -> list[Path]:
    """Load pygenomeviz example genbank dataset

    Load genbank datasets from <https://github.com/moshi4/pygenomeviz-data-v1>
    and cache datasets in local directory (Default: `~/.cache/pygenomeviz/`).

    List of dataset name

    - `acinetobacter_phage` (4 species)
    - `yersinia_phage` (4 species)
    - `enterobacteria_phage` (6 species)
    - `mycoplasma_mycoides` (4 species)
    - `escherichia_coli` (4 species, gzip compressed)
    - `saccharomyces` (3 species, gzip compressed)

    Parameters
    ----------
    name : str
        Dataset name (e.g. `enterobacteria_phage`)
    cache_dir : str | Path | None, optional
        Output cache directory (Default: `~/.cache/pygenomeviz/`)
    overwrite_cache : bool, optional
        If True, overwrite cached dataset

    Returns
    -------
    gbk_files : list[Path]
        Genbank files
    """
    logger = logging.getLogger(__name__)

    # Check specified name dataset exists or not
    if name not in GBK_DATASET.keys():
        raise ValueError(f"'{name}' dataset not found.")

    # Dataset cache local directory
    if cache_dir is None:
        package_name = __name__.split(".")[0]
        cache_base_dir = Path.home() / ".cache" / package_name
        cache_dir = cache_base_dir / "genbank" / name
        os.makedirs(cache_dir, exist_ok=True)
    else:
        cache_dir = Path(cache_dir)

    # Download & cache dataset
    gbk_files: list[Path] = []
    for filename in GBK_DATASET[name]:
        file_url = GITHUB_DATA_URL + f"genbank/{name}/{filename}"
        file_path = cache_dir / filename
        if overwrite_cache or not file_path.exists():
            logger.info(f"Download '{file_url}'...")
            urlretrieve(file_url, file_path)
            logger.info(f"Save download file in '{file_path}'")
        else:
            logger.info(f"Cached genbank file found in '{file_path}'")
        gbk_files.append(file_path)

    return gbk_files


def load_example_gff_file(
    filename: GffExampleFileName,
    *,
    cache_dir: str | Path | None = None,
    overwrite_cache: bool = False,
) -> Path:
    """Load pygenomeviz example GFF file

    Load example GFF file from <https://github.com/moshi4/pygenomeviz-data-v1/>
    and cache GFF file in local directory (Default: `~/.cache/pygenomeviz/`).

    List of example GFF filename

    - `enterobacteria_phage.gff`
    - `mycoplasma_mycoides.gff`
    - `escherichia_coli.gff.gz`
    - `saccharomyces_cerevisiae.gff.gz`

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
    if filename not in GFF_FILES:
        raise ValueError(f"{filename=} not found.")

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
    >>> gbk_text = fetch_genbank_by_accid("NC_013600")
    >>> gbk = Genbank(gbk_text)
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
        with open(gbk_outfile, "w", encoding="utf-8") as f:
            f.write(gbk_text)
        gbk_fetch_data = StringIO(gbk_text)  # type: ignore

    return gbk_fetch_data
