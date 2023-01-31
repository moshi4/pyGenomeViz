from __future__ import annotations

import argparse
import os
from pathlib import Path

from pygenomeviz import Genbank, GenomeViz, __version__
from pygenomeviz.align import AlignCoord, MMseqs
from pygenomeviz.config import LiteralTypes
from pygenomeviz.scripts import gbk_files2gbk_objects, get_argparser, print_args


def main():
    """Main function called from CLI"""
    # Get arguments
    args = get_args()
    print_args(args)
    # Run workflow
    run(**args.__dict__)


def run(
    # General options
    gbk_resources: list[str | Path] | list[Genbank],
    outdir: str | Path,
    format: list[str] = ["png", "html"],
    reuse: bool = False,
    # MMseqs options
    evalue: float = 1e-3,
    min_identity: float = 0,
    thread_num: int | None = None,
    # Figure appearence options
    fig_width: float = 15,
    fig_track_height: float = 1.0,
    feature_track_ratio: float = 1.0,
    link_track_ratio: float = 5.0,
    tick_track_ratio: float = 1.0,
    track_labelsize: int = 20,
    tick_labelsize: int = 15,
    normal_link_color: str = "grey",
    inverted_link_color: str = "red",
    align_type: LiteralTypes.ALIGN_TYPE = "center",
    tick_style: LiteralTypes.TICK_STYLE = None,
    feature_plotstyle: LiteralTypes.PLOTSTYLE = "bigarrow",
    arrow_shaft_ratio: float = 0.5,
    feature_color: str = "orange",
    feature_linewidth: float = 0,
    colorbar_width: float = 0.01,
    colorbar_height: float = 0.2,
    curve: bool = True,
    dpi: int = 300,
) -> GenomeViz:
    """Run genome visualization workflow using MMseqs

    Parameters
    ----------
    gbk_resources : list[str | Path] | list[Genbank]
        Input genome genbank files or Genbank objects
    outdir : str | Path
        Output directory
    format : list[str], optional
        Output image format (`png`|`jpg`|`svg`|`pdf`|`html`)
    reuse : bool, optional
        If True, reuse previous result if available
    evalue : float, optional
        MMseqs RBH search E-value parameter
    min_identity : float, optional
        Min-identity threshold to be plotted
    thread_num : int | None, optional
        MMseqs thread number to be used. If None, MaxThread - 1.
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
    feature_plotstyle : str, optional
        Feature plotstyle (`bigarrow`|`arrow`)
    arrow_shaft_ratio : float, optional
        Feature arrow shaft ratio
    feature_color : str, optional
        Feature color
    feature_linewidth : float, optional
        Feature edge line width
    colorbar_width : float, optional
        Colorbar width
    colorbar_height : float, optional
        Colorbar height
    curve : bool, optional
        If True, plot curved style link
    dpi : int, optional
        Figure DPI

    Returns
    -------
    gv : GenomeViz
        GenomeViz instance
    """
    # Check MMseqs installation
    MMseqs.check_installation()

    # Setup output contents
    outdir = Path(outdir)
    mmseqs_dir = outdir / "mmseqs_result"
    os.makedirs(mmseqs_dir, exist_ok=True)
    align_coords_file = outdir / "align_coords.tsv"

    # Setup Genbank objects
    gbk_list: list[Genbank] = []
    for gr in gbk_resources:
        if isinstance(gr, Genbank):
            gbk_list.append(gr)
        else:
            gbk_list.append(Genbank(gr))

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

    # Set tracks & features
    for gbk in gbk_list:
        track = gv.add_feature_track(
            gbk.name,
            size=gbk.range_size,
            start_pos=gbk.min_range,
            labelsize=track_labelsize,
        )
        track.add_genbank_features(
            gbk,
            feature_type="CDS",
            plotstyle=feature_plotstyle,
            arrow_shaft_ratio=arrow_shaft_ratio,
            facecolor=feature_color,
            linewidth=feature_linewidth,
            allow_partial=False,
        )

    # MMseqs alignment
    if align_coords_file.exists() and reuse:
        print("Reuse previous MMseqs RBH search result.\n")
        align_coords = AlignCoord.read(align_coords_file)
    else:
        print("Run MMseqs RBH search.\n")
        mmseqs = MMseqs(gbk_list, mmseqs_dir, 0, evalue, thread_num)
        align_coords = mmseqs.run()
        AlignCoord.write(align_coords, align_coords_file)
    align_coords = AlignCoord.filter(align_coords, 0, min_identity)

    # Set links
    min_identity = int(min([ac.identity for ac in align_coords]))
    for ac in align_coords:
        gv.add_link(
            ac.ref_link,
            ac.query_link,
            normal_link_color,
            inverted_link_color,
            curve=curve,
            v=ac.identity,
            vmin=min_identity,
        )

    # Plot figure
    fig = gv.plotfig(dpi=dpi)

    # Set colorbar
    bar_colors = [normal_link_color]
    contain_inverted_align = any([ac.is_inverted for ac in align_coords])
    if contain_inverted_align:
        bar_colors.append(inverted_link_color)
    gv.set_colorbar(
        fig,
        bar_colors,
        vmin=min_identity,
        bar_height=colorbar_height,
        bar_width=colorbar_width,
    )

    # Save figure
    for fmt in format:
        output_file = outdir / f"result.{fmt}"
        if fmt == "html":
            gv.savefig_html(output_file, fig)
        else:
            fig.savefig(output_file)
        print(f"Save {fmt} format result figure ({output_file}).")

    return gv


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
    parser = get_argparser(prog_name="MMseqs")

    #######################################################
    # General options
    #######################################################
    general_opts = parser.add_argument_group("General Options")
    general_opts.add_argument(
        "--gbk_resources",
        type=str,
        help="Input genome genbank file resources\n"
        "User can optionally specify genome range and reverse complement.\n"
        "- Example1. Set 100 - 1000 range 'file:100-1000'\n"
        "- Example2. Set reverse complement 'file::-1'\n"
        "- Example3. Set 100 - 1000 range of reverse complement 'file:100-1000:-1'",
        nargs="+",
        required=True,
        metavar="IN",
    )
    general_opts.add_argument(
        "-o",
        "--outdir",
        type=str,
        help="Output directory",
        required=True,
        metavar="OUT",
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
    # MMseqs options
    #######################################################
    mmseqs_opts = parser.add_argument_group("MMseqs Options")
    default_evalue = 1e-3
    mmseqs_opts.add_argument(
        "-e",
        "--evalue",
        type=float,
        help=f"MMseqs RBH search E-value parameter (Default: {default_evalue:.0e})",
        default=default_evalue,
        metavar="",
    )
    default_min_identity = 0
    mmseqs_opts.add_argument(
        "--min_identity",
        type=float,
        help=f"Min-identity threshold to be plotted (Default: {default_min_identity})",
        default=default_min_identity,
        metavar="",
    )
    cpu_num = os.cpu_count()
    default_thread_num = 1 if cpu_num is None or cpu_num == 1 else cpu_num - 1
    mmseqs_opts.add_argument(
        "-t",
        "--thread_num",
        type=int,
        help="Threads number parameter (Default: MaxThread - 1)",
        default=default_thread_num,
        metavar="",
    )
    #######################################################
    # Figure appearence options
    #######################################################
    fig_opts = parser.add_argument_group("Figure Appearence Options")
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
    default_inverted_link_color = "red"
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
    feature_plotstyle = "bigarrow"
    feature_plotstyle_list = ["bigarrow", "arrow"]
    fig_opts.add_argument(
        "--feature_plotstyle",
        type=str,
        help="Feature plot style ('bigarrow'[*]|'arrow')",
        default=feature_plotstyle,
        choices=feature_plotstyle_list,
        metavar="",
    )
    default_arrow_shaft_ratio = 0.5
    fig_opts.add_argument(
        "--arrow_shaft_ratio",
        type=float,
        help=f"Feature arrow shaft ratio (Default: {default_arrow_shaft_ratio})",
        default=default_arrow_shaft_ratio,
        metavar="",
    )
    default_feature_color = "orange"
    fig_opts.add_argument(
        "--feature_color",
        type=str,
        help=f"Feature color (Default: '{default_feature_color}')",
        default=default_feature_color,
        metavar="",
    )
    default_feature_linewidth = 0.0
    fig_opts.add_argument(
        "--feature_linewidth",
        type=float,
        help=f"Feature edge line width (Default: {default_feature_linewidth})",
        default=default_feature_linewidth,
        metavar="",
    )
    default_colorbar_width = 0.01
    fig_opts.add_argument(
        "--colorbar_width",
        type=float,
        help=f"Colorbar width (Default: {default_colorbar_width})",
        default=default_colorbar_width,
        metavar="",
    )
    default_colorbar_height = 0.2
    fig_opts.add_argument(
        "--colorbar_height",
        type=float,
        help=f"Colorbar height (Default: {default_colorbar_height})",
        default=default_colorbar_height,
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
    if len(args.gbk_resources) < 2:
        err_msg = "--gbk_resources must be set at least two files.\n"
        parser.error(err_msg)

    args.gbk_resources = gbk_files2gbk_objects(args.gbk_resources)

    return args


if __name__ == "__main__":
    main()
