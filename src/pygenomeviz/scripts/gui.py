#!/usr/bin/env python
from __future__ import annotations

import argparse
import importlib.util
import os
import subprocess as sp
import sys
import textwrap
from pathlib import Path

import pygenomeviz

CLI_NAME = "pgv-gui"


def main() -> None:
    """Launch pyGenomeViz WebApp"""
    # Get arguments
    args = get_args()
    port: int = args.port

    # Check streamlit installation
    if importlib.util.find_spec("streamlit") is None:
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
        sys.exit(1)

    # Launch app
    app_path = Path(__file__).parent.parent / "gui" / "app.py"
    os.environ["STREAMLIT_THEME_BASE"] = "dark"
    os.environ["STREAMLIT_SERVER_RUN_ON_SAVE"] = "true"
    os.environ["STREAMLIT_BROWSER_GATHER_USAGE_STATS"] = "false"
    os.environ["STREAMLIT_SERVER_MAX_UPLOAD_SIZE"] = "100"
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
    description = textwrap.dedent(
        """
        pyGenomeViz CLI for launching Streamlit Web Application

        Users can access the web app with http://localhost:8501 (default).
        """
    )[1:-1]
    parser = argparse.ArgumentParser(
        description=description,
        usage=f"{CLI_NAME} [options]",
        epilog="",
        add_help=False,
        formatter_class=argparse.RawTextHelpFormatter,
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

    args = parser.parse_args(cli_args)

    return args
