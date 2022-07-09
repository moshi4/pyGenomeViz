import argparse
import os
import shutil
from pathlib import Path

from pygenomeviz import __version__
from pygenomeviz.utils import DATASETS, load_dataset


def main():
    """Download pyGenomeViz dataset"""
    # Get arguments
    args = get_args()
    dataset_name: str = args.name
    outdir: Path = args.outdir

    # Download dataset
    os.makedirs(outdir, exist_ok=True)
    gbk_files, _ = load_dataset(dataset_name)
    for gbk_file in gbk_files:
        shutil.copy(gbk_file, outdir)


def get_args() -> argparse.Namespace:
    """Get arguments

    Returns
    -------
    args : argparse.Namespace
        Argument parameters
    """
    description = "Download pyGenomeViz genbank dataset"
    dataset_names = list(DATASETS.keys())
    epilog = f"Dataset name = {dataset_names}"
    parser = argparse.ArgumentParser(
        description=description, epilog=epilog, add_help=False
    )

    parser.add_argument(
        "-n",
        "--name",
        type=str,
        help="Dataset name",
        required=True,
        choices=dataset_names,
        metavar="",
    )
    default_outdir = Path("./")
    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        help=f"Output directory (Default: {default_outdir})",
        default=default_outdir,
        metavar="",
    )
    parser.add_argument(
        "-v",
        "--version",
        version=f"v{__version__}",
        help="Print version information",
        action="version",
    )
    parser.add_argument(
        "-h",
        "--help",
        help="Show this help message and exit",
        action="help",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
