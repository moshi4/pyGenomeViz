#!/usr/bin/env python
from __future__ import annotations

import argparse
import os
import shutil
from pathlib import Path

import pygenomeviz
from pygenomeviz.logger import get_logger
from pygenomeviz.scripts import CustomHelpFormatter
from pygenomeviz.typing import GenbankDatasetName
from pygenomeviz.utils import load_example_genbank_dataset
from pygenomeviz.utils.download import GBK_DATASET

CLI_NAME = "pgv-download"


def main():
    """Main function called from CLI"""
    # Get arguments
    args = get_args()
    dataset_name: GenbankDatasetName = args.dataset_name
    outdir: Path = args.outdir
    quiet: bool = args.quiet
    cache_only: bool = args.cache_only

    logger = get_logger(__name__, quiet=quiet)

    # Download dataset
    if cache_only:
        load_example_genbank_dataset(dataset_name, quiet=quiet)
    else:
        os.makedirs(outdir, exist_ok=True)
        gbk_files = load_example_genbank_dataset(dataset_name, quiet=quiet)
        for gbk_file in gbk_files:
            logger.info(f"Copy '{gbk_file}' to '{outdir}/{gbk_file.name}'")
            shutil.copy(gbk_file, outdir)


def get_args() -> argparse.Namespace:
    """Get arguments

    Returns
    -------
    args : argparse.Namespace
        Argument parameters
    """
    dataset_names = tuple(GBK_DATASET.keys())
    parser = argparse.ArgumentParser(
        description="pyGenomeViz CLI for downloading example genbank dataset",
        usage=f"{CLI_NAME} [dataset_name] -o [outdir]",
        epilog=f"dataset_name = {dataset_names}",
        add_help=False,
        allow_abbrev=False,
        formatter_class=CustomHelpFormatter,
    )

    parser.add_argument(
        "dataset_name",
        nargs=1,
        type=str,
        help="Target dataset name",
    )

    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        help="Output directory (Default: '.')",
        default=Path("./"),
        metavar="",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        help="No print log on screen (default: OFF)",
        action="store_true",
    )
    parser.add_argument(
        "--cache_only",
        help=argparse.SUPPRESS,
        action="store_true",
    )
    parser.add_argument(
        "-v",
        "--version",
        version=f"v{pygenomeviz.__version__}",
        help="Print version information",
        action="version",
    )
    parser.add_argument(
        "-h",
        "--help",
        help="Show this help message and exit",
        action="help",
    )

    # Validate dataset name argument
    args = parser.parse_args()
    dataset_name = str(args.dataset_name[0])
    if dataset_name not in dataset_names:
        parser.error(f"{dataset_name=} is invalid.\nChoose from {dataset_names}.")
    args.dataset_name = dataset_name

    return args


if __name__ == "__main__":
    main()
