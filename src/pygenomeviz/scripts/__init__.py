from __future__ import annotations

import argparse
import logging
import os
import platform
import sys
from pathlib import Path

import Bio
import matplotlib
from matplotlib.colors import is_color_like

import pygenomeviz
from pygenomeviz.align import Blast, MMseqs, MUMmer, ProgressiveMauve
from pygenomeviz.typing import AlnCliName

LOG_FILENAME = "pgv-cli.log"
ALIGN_COORDS_FILENAME = "align_coords.tsv"

CLI_NAME2TOOL_NAME: dict[AlnCliName, str] = {
    "pgv-blast": Blast.get_tool_name(),
    "pgv-mummer": MUMmer.get_tool_name(),
    "pgv-mmseqs": MMseqs.get_tool_name(),
    "pgv-pmauve": ProgressiveMauve.get_tool_name(),
}


class CustomHelpFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, prog, indent_increment=2, max_help_position=40, width=None):
        super().__init__(prog, indent_increment, max_help_position, width)

    def _format_args(self, action, default_metavar):
        return ""


def log_basic_env_info(
    logger: logging.Logger,
    cli_name: AlnCliName,
    log_params: dict | None,
) -> None:
    """Logging basic environment information

    Parameters
    ----------
    logger : logging.Logger
        Logger
    cli_name : str
        CLI name
    log_params : dict | None
        Log parameters
    """
    logger.info(f"Run pyGenomeViz v{pygenomeviz.__version__} CLI workflow ({cli_name})")
    logger.info(f"$ {Path(sys.argv[0]).name} {' '.join(sys.argv[1:])}")
    logger.info(f"Operating System: {sys.platform}")
    logger.info(f"Python Version: v{platform.python_version()}")
    logger.info(f"Check Dependencies: matplotlib v{matplotlib.__version__}")  # type: ignore
    logger.info(f"Check Dependencies: biopython v{Bio.__version__}")
    if log_params:
        for k, v in log_params.items():
            logger.info(f"Parameters: {k}={repr(v)}")


def setup_argparser(
    parser: argparse.ArgumentParser,
    cli_name: AlnCliName,
) -> None:
    """Setup argument parser

    Parameters
    ----------
    parser : argparse.ArgumentParser
        Argument parser
    cli_name : CliNameList
        Setup target command line name
    """
    # Positional arguments
    if cli_name == "pgv-pmauve":
        seqs_help_msg = "Input genbank or fasta files"
    else:
        seqs_help_msg = "Input genbank files"

    parser.add_argument(
        "seqs",
        nargs="+",
        type=str,
        help=seqs_help_msg,
    )

    # General Options
    general_arg_group = parser.add_argument_group("General Options")
    _setup_general_arg_group(general_arg_group)

    # Alignment Options
    if cli_name in ("pgv-blast", "pgv-mummer", "pgv-mmseqs"):
        tool_name = CLI_NAME2TOOL_NAME[cli_name]
        align_arg_group = parser.add_argument_group(f"{tool_name} Alignment Options")
        _setup_align_arg_group(align_arg_group, cli_name)

    # Figure Appearence Options
    fig_arg_group = parser.add_argument_group("Figure Appearence Options")
    _setup_fig_arg_group(fig_arg_group, cli_name)


def _setup_general_arg_group(general_arg_group: argparse._ArgumentGroup) -> None:
    """Setup general argument group for pgv-[blast|mummer|mmseqs|pmauve]

    Parameters
    ----------
    general_arg_group : argparse._ArgumentGroup
        General argument group
    """
    general_arg_group.add_argument(
        "-o",
        "--outdir",
        type=str,
        help="Output directory",
        required=True,
        metavar="",
    )
    default_formats = ("png", "html")
    general_arg_group.add_argument(
        "--formats",
        type=str,
        nargs="+",
        help="Output image format ('png'[*],'jpg','svg','pdf',`html`[*])",
        default=default_formats,
        choices=("png", "jpg", "svg", "pdf", "html"),
        metavar="",
    )
    general_arg_group.add_argument(
        "--reuse",
        help="Reuse previous alignment result if available",
        action="store_true",
    )
    general_arg_group.add_argument(
        "-q",
        "--quiet",
        help="No print log on screen (default: OFF)",
        action="store_true",
    )
    general_arg_group.add_argument(
        "--debug",
        help=argparse.SUPPRESS,
        action="store_true",
    )
    general_arg_group.add_argument(
        "-v",
        "--version",
        version=f"v{pygenomeviz.__version__}",
        help="Print version information",
        action="version",
    )
    general_arg_group.add_argument(
        "-h",
        "--help",
        help="Show this help message and exit",
        action="help",
    )


def _setup_align_arg_group(
    align_arg_group: argparse._ArgumentGroup,
    cli_name: AlnCliName,
) -> None:
    """Setup alignment argument group for pgv-[blast|mummer|mmseqs]

    Parameters
    ----------
    align_arg_group : argparse._ArgumentGroup
        Alignment argument group
    cli_name : CliName
        Target command line name
    """
    if cli_name in ("pgv-blast", "pgv-mummer"):
        default_seqtype = "nucleotide"
        align_arg_group.add_argument(
            "--seqtype",
            type=str,
            help="Alignment sequence type ('nucleotide'[*]|'protein')",
            default=default_seqtype,
            choices=("nucleotide", "protein"),
            metavar="",
        )
    cpu_num = os.cpu_count()
    default_threads = 1 if cpu_num is None or cpu_num == 1 else cpu_num - 1
    align_arg_group.add_argument(
        "--threads",
        type=int,
        help="Threads number (Default: MaxThread - 1)",
        default=default_threads,
        metavar="",
    )
    default_length_thr = 0
    align_arg_group.add_argument(
        "--length_thr",
        type=int,
        help=f"Length threshold to be plotted (Default: {default_length_thr})",
        default=default_length_thr,
        metavar="",
    )
    default_identity_thr = 0
    align_arg_group.add_argument(
        "--identity_thr",
        type=float,
        help=f"Identity threshold to be plotted (Default: {default_identity_thr})",
        default=default_identity_thr,
        metavar="",
    )
    if cli_name in ("pgv-blast", "pgv-mmseqs"):
        default_evalue_thr = 1e-3
        align_arg_group.add_argument(
            "--evalue_thr",
            type=float,
            help=f"E-value threshold to be plotted (Default: {default_evalue_thr:.0e})",
            default=default_evalue_thr,
            metavar="",
        )


def _setup_fig_arg_group(
    fig_arg_group: argparse._ArgumentGroup,
    cli_name: AlnCliName,
) -> None:
    """Setup figure appearence argument group for pgv-[blast|mummer|mmseqs|pmauve]

    Parameters
    ----------
    fig_arg_group : argparse._ArgumentGroup
        Figure appearence argument group
    cli_name : CliName
        Target command line name
    """
    default_fig_width = 15
    fig_arg_group.add_argument(
        "--fig_width",
        type=float,
        help=f"Figure width (Default: {default_fig_width})",
        default=default_fig_width,
        metavar="",
    )
    default_fig_track_height = 1.0
    fig_arg_group.add_argument(
        "--fig_track_height",
        type=float,
        help=f"Figure track height (Default: {default_fig_track_height})",
        default=default_fig_track_height,
        metavar="",
    )
    default_track_align_type = "center"
    fig_arg_group.add_argument(
        "--track_align_type",
        type=str,
        help="Figure tracks align type ('left'|'center'[*]|'right')",
        default=default_track_align_type,
        choices=("left", "center", "right"),
        metavar="",
    )
    default_feature_track_ratio = 0.25
    fig_arg_group.add_argument(
        "--feature_track_ratio",
        type=float,
        help=f"Feature track ratio (Default: {default_feature_track_ratio})",
        default=default_feature_track_ratio,
        metavar="",
    )
    fig_arg_group.add_argument(
        "--show_scale_bar",
        help="Show scale bar (Default: OFF)",
        action="store_true",
    )
    fig_arg_group.add_argument(
        "--show_scale_xticks",
        help="Show scale xticks (Default: OFF)",
        action="store_true",
    )
    fig_arg_group.add_argument(
        "--curve",
        help="Plot curved style link (Default: OFF)",
        action="store_true",
    )
    default_dpi = 300
    fig_arg_group.add_argument(
        "--dpi",
        type=int,
        help=f"Figure DPI (Default: {default_dpi})",
        default=default_dpi,
        metavar="",
    )
    default_track_labelsize = 20
    fig_arg_group.add_argument(
        "--track_labelsize",
        type=int,
        help=f"Track label size (Default: {default_track_labelsize})",
        default=default_track_labelsize,
        metavar="",
    )
    default_scale_labelsize = 15
    fig_arg_group.add_argument(
        "--scale_labelsize",
        type=int,
        help=f"Scale label size (Default: {default_scale_labelsize})",
        default=default_scale_labelsize,
        metavar="",
    )
    default_normal_link_color = "grey"
    fig_arg_group.add_argument(
        "--normal_link_color",
        type=str,
        help=f"Normal link color (Default: '{default_normal_link_color}')",
        default=default_normal_link_color,
        metavar="",
    )
    default_inverted_link_color = "red"
    fig_arg_group.add_argument(
        "--inverted_link_color",
        type=str,
        help=f"Inverted link color (Default: '{default_inverted_link_color}')",
        default=default_inverted_link_color,
        metavar="",
    )

    if cli_name in ("pgv-blast", "pgv-mummer", "pgv-mmseqs"):
        default_segment_space = 0.02
        fig_arg_group.add_argument(
            "--segment_space",
            type=float,
            help=f"Track segment space ratio (Default: {default_segment_space})",
            default=default_segment_space,
            metavar="",
        )
        default_feature_type2color = ["CDS:orange"]
        fig_arg_group.add_argument(
            "--feature_type2color",
            nargs="*",
            type=str,
            help=f"Feature plot type & color (Default: {default_feature_type2color})",
            default=default_feature_type2color,
            metavar="",
        )
        default_pseudo_color = "lightgrey"
        fig_arg_group.add_argument(
            "--pseudo_color",
            type=str,
            help=f"Pseudo feature plot color (Default: '{default_pseudo_color}')",
            default=default_pseudo_color,
            metavar="",
        )
        default_feature_plotstyle = "arrow"
        fig_arg_group.add_argument(
            "--feature_plotstyle",
            type=str,
            help="Feature plot style ('[big]arrow'[*]|'[big]box'|'[big]rbox')",
            default=default_feature_plotstyle,
            choices=("bigarrow", "arrow", "bigbox", "box", "bigrbox", "rbox"),
            metavar="",
        )
        default_feature_linewidth = 0.0
        fig_arg_group.add_argument(
            "--feature_linewidth",
            type=float,
            help=f"Feature line width (Default: {default_feature_linewidth})",
            default=default_feature_linewidth,
            metavar="",
        )
        default_feature_labeltrack = "top"
        fig_arg_group.add_argument(
            "--feature_labeltrack",
            type=str,
            help="Feature label target track ('top'[*]|'all')",
            default=default_feature_labeltrack,
            choices=("top", "all"),
            metavar="",
        )
        default_feature_labeltype = "None"
        fig_arg_group.add_argument(
            "--feature_labeltype",
            type=lambda v: None if v == "None" else v,
            help="Feature label type ('product'|'gene'|'protein_id'|'None'[*])",
            default=default_feature_labeltype,
            choices=("product", "gene", "protein_id", None),
            metavar="",
        )
        default_feature_labelsize = 8
        fig_arg_group.add_argument(
            "--feature_labelsize",
            type=int,
            help=f"Feature label size (Default: {default_feature_labelsize})",
            default=default_feature_labelsize,
            metavar="",
        )
        default_cbar_width = 0.01
        fig_arg_group.add_argument(
            "--cbar_width",
            type=float,
            help=f"Colorbar width (Default: {default_cbar_width})",
            default=default_cbar_width,
            metavar="",
        )
        default_cbar_height = 0.2
        fig_arg_group.add_argument(
            "--cbar_height",
            type=float,
            help=f"Colorbar height (Default: {default_cbar_height})",
            default=default_cbar_height,
            metavar="",
        )

    if cli_name == "pgv-pmauve":
        default_refid = 0
        fig_arg_group.add_argument(
            "--refid",
            type=int,
            help=f"Reference genome index (Default: {default_refid})",
            default=default_refid,
            metavar="",
        )
        default_block_plotstyle = "box"
        fig_arg_group.add_argument(
            "--block_plotstyle",
            type=str,
            help="Synteny block plot style ('box'[*]|'bigbox')",
            default=default_block_plotstyle,
            choices=("box", "bigbox"),
            metavar="",
        )
        default_block_cmap = "gist_rainbow"
        fig_arg_group.add_argument(
            "--block_cmap",
            type=str,
            help=f"Synteny block colormap (Default: '{default_block_cmap}')",
            default=default_block_cmap,
            metavar="",
        )


def validate_args(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    """Validate command arguments

    Parameters
    ----------
    args : argparse.Namespace
        Command arguments
    parser : argparse.ArgumentParser
        Command parser
    """
    # seqs (target = all)
    if hasattr(args, "seqs"):
        if len(args.seqs) < 2:
            parser.error("Input must be set at least two files.")

    # normal_link_color (target = all)
    if hasattr(args, "normal_link_color"):
        if not is_color_like(args.normal_link_color):
            parser.error(f"{args.normal_link_color=} is invalid color name or code.")
    # inverted_link_color (target = all)
    if hasattr(args, "inverted_link_color"):
        if not is_color_like(args.inverted_link_color):
            parser.error(f"{args.inverted_link_color=} is invalid color name or code.")

    # feature_type2color (target = blast, mummer, mmseqs)
    if hasattr(args, "feature_type2color"):
        feature_type2color = {}
        for feature_type_and_color in args.feature_type2color:
            try:
                feature_type, color = str(feature_type_and_color).split(":")
            except Exception:
                parser.error(f"Failed to parse args ({args.feature_type2color=}).")
            if not is_color_like(color):
                parser.error(f"{color=} is invalid ({args.feature_type2color=}).")
            feature_type2color[feature_type] = color
        args.feature_type2color = feature_type2color

    # pseudo_color (target = blast, mummer, mmseqs)
    if hasattr(args, "pseudo_color"):
        if not is_color_like(args.pseudo_color):
            parser.error(f"{args.pseudo_color=} is invalid color name or code.")

    # refid (target = pmauve)
    if hasattr(args, "refid"):
        min_refid, max_refid = 0, len(args.seqs) - 1
        if not min_refid <= args.refid <= max_refid:
            parser.error(f"--refid must be '{min_refid} <= refid <= {max_refid}'.")

    # cmap (target = pmauve)
    if hasattr(args, "cmap"):
        colormaps = list(matplotlib.colormaps)  # type: ignore
        if args.cmap not in colormaps:
            parser.error(f"{args.cmap=} is invalid colormap.\nAvailable {colormaps=}")
