from __future__ import annotations

import argparse
import os
import subprocess as sp
import textwrap
from importlib.util import find_spec
from pathlib import Path

from pygenomeviz import __version__


def main() -> None:
    """Launch pyGenomeViz WebApp"""
    # Get arguments
    args = get_args()
    port: int = args.port

    # Check streamlit installation
    if find_spec("streamlit") is None:
        err_msg = textwrap.dedent(
            """
            Failed to launch pyGenomeViz WebApp. Streamlit is not installed!!

            # Install PyPI package
            $ pip install -U pygenomeviz[gui]

            # Install conda-forge package
            $ conda install -c conda-forge streamlit
            """
        )
        print(err_msg)
        exit(1)

    # Launch app
    app_path = Path(__file__).parent.parent / "gui" / "app.py"
    os.environ["STREAMLIT_THEME_BASE"] = "dark"
    os.environ["STREAMLIT_SERVER_RUN_ON_SAVE"] = "true"
    os.environ["STREAMLIT_BROWSER_GATHER_USAGE_STATS"] = "false"
    os.environ["PGV_GUI_LOCAL"] = "true"
    sp.run(f"streamlit run {app_path} --server.port {port}", shell=True)


def get_args(cli_args: list[str] | None = None) -> argparse.Namespace:
    """Get arguments

    Parameters
    ----------
    cli_args : list[str] | None, optional
        CLI arguments (Used in unittest)

    Returns
    -------
    args : argparse.Namespace
        Argument parameters
    """
    description = "Launch pyGenomeViz WebApp"
    epilog = ""
    parser = argparse.ArgumentParser(
        description=description, epilog=epilog, add_help=False
    )

    default_port = 8501
    parser.add_argument(
        "-p",
        "--port",
        type=int,
        help=f"Port number to open web browser (Default: {default_port})",
        default=default_port,
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

    args = parser.parse_args(cli_args)

    return args
