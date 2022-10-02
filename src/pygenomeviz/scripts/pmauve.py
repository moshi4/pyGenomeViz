from __future__ import annotations

import argparse
import csv
import os
from collections import defaultdict
from pathlib import Path

from pygenomeviz import GenomeViz, __version__
from pygenomeviz.align import AlignCoord, ProgressiveMauve
from pygenomeviz.config import LiteralTypes
from pygenomeviz.scripts import get_argparser, print_args
from pygenomeviz.utils import ColorCycler


def main():
    """Main function called from CLI"""
    # Get arguments
    args = get_args()
    print_args(args)
    # Run workflow
    run(**args.__dict__)


def run(
    # General options
    seq_files: list[str | Path],
    outdir: str | Path,
    refid: int = 0,
    format: list[str] = ["png", "html"],
    reuse: bool = False,
    # Figure appearence options
    fig_width: float = 15,
    fig_track_height: float = 1.0,
    feature_track_ratio: float = 1.0,
    link_track_ratio: float = 5.0,
    tick_track_ratio: float = 1.0,
    track_labelsize: int = 20,
    tick_labelsize: int = 15,
    normal_link_color: str = "grey",
    inverted_link_color: str = "tomato",
    align_type: LiteralTypes.ALIGN_TYPE = "center",
    tick_style: LiteralTypes.TICK_STYLE = None,
    plotstyle: LiteralTypes.PLOTSTYLE = "box",
    cmap: str = "hsv",
    curve: bool = True,
    dpi: int = 300,
) -> GenomeViz:
    """Run genome visualization workflow using progressiveMauve

    Parameters
    ----------
    seq_files : list[str | Path]
        Input genome sequence files (Genbank or Fasta format)
    outdir : str | Path
        Output directory
    format : list[str], optional
        Output image format (`png`|`jpg`|`svg`|`pdf`|`html`)
    reuse : bool, optional
        If True, reuse previous result if available
    fig_width : float, optional
        Figure width
    fig_track_height : float, optional
        Figure track height
    feature_track_ratio : float, optional
        Feature track ratio
    link_track_ratio : float, optional
        Link track ratio
    tick_track_ratio : float, optional
        Tick track ratio
    track_labelsize : int, optional
        Track label size
    tick_labelsize : int, optional
        Tick label size
    normal_link_color : str, optional
        Normal link color
    inverted_link_color : str, optional
        Inverted link color
    align_type : str, optional
        Figure tracks align type (`left`|`center`|`right`)
    tick_style : str | None, optional
        Tick style (`bar`|`axis`|`None`)
    plotstyle : str, optional
        Block box plot style (`box`|`bigbox`)
    cmap : str, optional
        Block box colormap
    curve : bool, optional
        If True, plot curved style link
    dpi : int, optional
        Figure DPI

    Returns
    -------
    gv : GenomeViz
        GenomeViz instance
    """
    # Check progressiveMauve installation
    ProgressiveMauve.check_installation()

    # Setup output contents
    outdir = Path(outdir)
    pmauve_dir = outdir / "pmauve_result"
    os.makedirs(pmauve_dir, exist_ok=True)
    align_coords_file = outdir / "align_coords.tsv"
    ColorCycler.set_cmap(cmap)

    # progressiveMauve alignment
    pmauve = ProgressiveMauve(seq_files, pmauve_dir, refid)
    if pmauve.bbone_file.exists() and reuse:
        print("Reuse previous progressiveMauve result.\n")
        align_coords = pmauve.parse_pmauve_file(pmauve.bbone_file)
        AlignCoord.write(align_coords, align_coords_file)
    else:
        print("Run progressiveMauve alignment.\n")
        align_coords = pmauve.run()
        AlignCoord.write(align_coords, align_coords_file)
    name2maxsize = get_name2maxsize(pmauve.filenames, pmauve.bbone_file)
    name2blocks = get_name2blocks(align_coords)

    # Setup GenomeViz instance
    gv = GenomeViz(
        fig_width=fig_width,
        fig_track_height=fig_track_height,
        feature_track_ratio=feature_track_ratio,
        link_track_ratio=link_track_ratio,
        tick_track_ratio=tick_track_ratio,
        align_type=align_type,
        tick_style=tick_style,
        tick_labelsize=tick_labelsize,
        plot_size_thr=0,
    )

    # Set tracks & block features
    for name in pmauve.filenames:
        maxsize = name2maxsize[name]
        track = gv.add_feature_track(name, maxsize, labelsize=track_labelsize)
        blocks = name2blocks[name]
        colors = ColorCycler.get_color_list(len(blocks))
        for block, color in zip(blocks, colors):
            start, end, strand = block
            track.add_feature(start, end, strand, plotstyle=plotstyle, facecolor=color)

    # Set links
    for ac in align_coords:
        gv.add_link(
            ac.ref_link,
            ac.query_link,
            normal_link_color,
            inverted_link_color,
            curve=curve,
        )

    # Save figure
    for fmt in format:
        output_file = outdir / f"result.{fmt}"
        if fmt == "html":
            gv.savefig_html(output_file)
        else:
            gv.savefig(output_file, dpi=dpi)
        print(f"Save {fmt} format result figure ({output_file}).")

    return gv


def get_name2maxsize(names: list[str], bbone_file: str | Path) -> dict[str, int]:
    """Get name to maxsize genome dict from bbone file

    Parameters
    ----------
    names: list[str]
        Name labels
    bbone_file : str | Path
        progressiveMauve bbone format file

    Returns
    -------
    name2maxsize : dict[str, int]
        name to maxsize dict
    """
    name2maxsize = defaultdict(int)
    with open(bbone_file) as f:
        reader = csv.reader(f, delimiter="\t")
        header_row = next(reader)
        genome_num = int(len(header_row) / 2)
        for row in reader:
            row = [abs(int(col)) for col in row]
            for i in range(genome_num):
                name, size = names[i], max(row[i * 2], row[i * 2 + 1])
                if name2maxsize[name] < size:
                    name2maxsize[name] = size
    return name2maxsize


def get_name2blocks(
    align_coords: list[AlignCoord],
) -> dict[str, list[tuple[int, int, int]]]:
    """Get name to align coord blocks from align coords

    Parameters
    ----------
    align_coords : list[AlignCoord]
        Alignment coords

    Returns
    -------
    name2blocks : dict[str, list[tuple[int, int, int]]]
        name to align coord blocks
    """
    name2blocks = defaultdict(list)
    for ac in align_coords:
        if ac.ref_block not in name2blocks[ac.ref_name]:
            name2blocks[ac.ref_name].append(ac.ref_block)
        if ac.query_block not in name2blocks[ac.query_name]:
            name2blocks[ac.query_name].append(ac.query_block)
    return name2blocks


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
    parser = get_argparser(prog_name="progressiveMauve")

    #######################################################
    # General options
    #######################################################
    general_opts = parser.add_argument_group("General Options")
    general_opts.add_argument(
        "--seq_files",
        type=Path,
        help="Input genome sequence files (Genbank or Fasta format)",
        nargs="+",
        required=True,
        metavar="IN",
    )
    general_opts.add_argument(
        "-o",
        "--outdir",
        type=Path,
        help="Output directory",
        required=True,
        metavar="OUT",
    )
    default_refid = 0
    general_opts.add_argument(
        "--refid",
        type=int,
        help=f"Reference genome index (Default: {default_refid})",
        default=default_refid,
        metavar="",
    )
    default_format = ["png", "html"]
    format_list = ["png", "jpg", "svg", "pdf", "html"]
    general_opts.add_argument(
        "--format",
        type=str,
        nargs="+",
        help="Output image format ('png'[*]|'jpg'|'svg'|'pdf'|`html`[*])",
        default=default_format,
        choices=format_list,
        metavar="",
    )
    general_opts.add_argument(
        "--reuse",
        help="Reuse previous result if available",
        action="store_true",
    )
    general_opts.add_argument(
        "-v",
        "--version",
        version=f"v{__version__}",
        help="Print version information",
        action="version",
    )
    general_opts.add_argument(
        "-h",
        "--help",
        help="Show this help message and exit",
        action="help",
    )
    #######################################################
    # Figure appearence options
    #######################################################
    fig_opts = parser.add_argument_group("Figure Appearence options")
    default_fig_width = 15
    fig_opts.add_argument(
        "--fig_width",
        type=float,
        help=f"Figure width (Default: {default_fig_width})",
        default=default_fig_width,
        metavar="",
    )
    default_fig_track_height = 1.0
    fig_opts.add_argument(
        "--fig_track_height",
        type=float,
        help=f"Figure track height (Default: {default_fig_track_height})",
        default=default_fig_track_height,
        metavar="",
    )
    default_feature_track_ratio = 1.0
    fig_opts.add_argument(
        "--feature_track_ratio",
        type=float,
        help=f"Feature track ratio (Default: {default_feature_track_ratio})",
        default=default_feature_track_ratio,
        metavar="",
    )
    default_link_track_ratio = 5.0
    fig_opts.add_argument(
        "--link_track_ratio",
        type=float,
        help=f"Link track ratio (Default: {default_link_track_ratio})",
        default=default_link_track_ratio,
        metavar="",
    )
    default_tick_track_ratio = 1.0
    fig_opts.add_argument(
        "--tick_track_ratio",
        type=float,
        help=f"Tick track ratio (Default: {default_tick_track_ratio})",
        default=default_tick_track_ratio,
        metavar="",
    )
    default_track_labelsize = 20
    fig_opts.add_argument(
        "--track_labelsize",
        type=int,
        help=f"Track label size (Default: {default_track_labelsize})",
        default=default_track_labelsize,
        metavar="",
    )
    default_tick_labelsize = 15
    fig_opts.add_argument(
        "--tick_labelsize",
        type=int,
        help=f"Tick label size (Default: {default_tick_labelsize})",
        default=default_tick_labelsize,
        metavar="",
    )
    default_normal_link_color = "grey"
    fig_opts.add_argument(
        "--normal_link_color",
        type=str,
        help=f"Normal link color (Default: '{default_normal_link_color}')",
        default=default_normal_link_color,
        metavar="",
    )
    default_inverted_link_color = "tomato"
    fig_opts.add_argument(
        "--inverted_link_color",
        type=str,
        help=f"Inverted link color (Default: '{default_inverted_link_color}')",
        default=default_inverted_link_color,
        metavar="",
    )
    default_align_type = "center"
    align_type_list = ["left", "center", "right"]
    fig_opts.add_argument(
        "--align_type",
        type=str,
        help="Figure tracks align type ('left'|'center'[*]|'right')",
        default=default_align_type,
        choices=align_type_list,
        metavar="",
    )
    default_tick_style = None
    fig_opts.add_argument(
        "--tick_style",
        type=str,
        help="Tick style ('bar'|'axis'|None[*])",
        default=default_tick_style,
        choices=["bar", "axis"],
        metavar="",
    )
    plotstyle = "box"
    plotstyle_list = ["box", "bigbox"]
    fig_opts.add_argument(
        "--plotstyle",
        type=str,
        help="Block box plot style ('box'[*]|'bigbox')",
        default=plotstyle,
        choices=plotstyle_list,
        metavar="",
    )
    default_cmap = "hsv"
    fig_opts.add_argument(
        "--cmap",
        type=str,
        help=f"Block box colormap (Default: '{default_cmap}')",
        default=default_cmap,
        metavar="",
    )
    fig_opts.add_argument(
        "--curve",
        help="Plot curved style link (Default: OFF)",
        action="store_true",
    )
    default_dpi = 300
    fig_opts.add_argument(
        "--dpi",
        type=int,
        help=f"Figure DPI (Default: {default_dpi})",
        default=default_dpi,
        metavar="",
    )
    args = parser.parse_args(cli_args)

    # Check arguments
    err_info = ""
    if len(args.seq_files) < 2:
        err_info += "--seq_files must be set at least two files.\n"
    for file in args.seq_files:
        if not file.exists():
            err_info += f"File not found '{file}'.\n"
    min_refid, max_refid = 0, len(args.seq_files) - 1
    if not min_refid <= args.refid <= max_refid:
        err_info += f"--refid must be '{min_refid} <= refid <= {max_refid}'\n"

    if err_info != "":
        parser.error("\n" + err_info)

    return args


if __name__ == "__main__":
    main()
