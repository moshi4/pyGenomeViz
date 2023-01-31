from __future__ import annotations

import argparse
from pathlib import Path

from pygenomeviz import Genbank, GenomeViz, __version__
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
    outfile: str | Path,
    # Figure appearence options
    fig_width: float = 15,
    fig_track_height: float = 1.0,
    feature_track_ratio: float = 1.0,
    space_track_ratio: float = 5.0,
    tick_track_ratio: float = 1.0,
    track_labelsize: int = 20,
    tick_labelsize: int = 15,
    align_type: LiteralTypes.ALIGN_TYPE = "left",
    tick_style: LiteralTypes.TICK_STYLE = None,
    feature_plotstyle: LiteralTypes.PLOTSTYLE = "bigarrow",
    arrow_shaft_ratio: float = 0.5,
    feature_color: str = "orange",
    feature_linewidth: float = 0,
    dpi: int = 300,
) -> GenomeViz:
    """Run simple genome visualization workflow

    Parameters
    ----------
    gbk_resources : list[str | Path] | list[Genbank]
        Input genome genbank files or Genbank objects
    outfile : str | Path
        Output file ('*.png'|'*.jpg'|'*.svg'|'*.pdf'|'*.html')
    fig_width : float, optional
        Figure width
    fig_track_height : float, optional
        Figure track height
    feature_track_ratio : float, optional
        Feature track ratio
    space_track_ratio : float, optional
        Space track ratio (track between feature tracks)
    tick_track_ratio : float, optional
        Tick track ratio
    track_labelsize : int, optional
        Track label size
    tick_labelsize : int, optional
        Tick label size
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
    dpi : int, optional
        Figure DPI

    Returns
    -------
    gv : GenomeViz
        GenomeViz instance
    """
    outfile = Path(outfile)

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
        link_track_ratio=space_track_ratio,
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

    # Save figure
    file_format = outfile.suffix.replace(".", "")
    if file_format == "html":
        gv.savefig_html(outfile)
    else:
        gv.savefig(outfile, dpi=dpi)
    print(f"Save {file_format} format result figure ({outfile}).")

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
    parser = get_argparser()

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
        "--outfile",
        type=str,
        help="Output file ('*.png'|'*.jpg'|'*.svg'|'*.pdf'|'*.html')",
        required=True,
        metavar="OUT",
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
    default_space_track_ratio = 2.0
    fig_opts.add_argument(
        "--space_track_ratio",
        type=float,
        help=f"Space (b/w feature) track ratio (Default: {default_space_track_ratio})",
        default=default_space_track_ratio,
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
    default_align_type = "left"
    align_type_list = ["left", "center", "right"]
    fig_opts.add_argument(
        "--align_type",
        type=str,
        help="Figure tracks align type ('left'[*]|'center'|'right')",
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
    format_list = [".png", ".jpg", ".svg", ".pdf", ".html"]
    if Path(args.outfile).suffix not in format_list:
        err_msg = f"--outfile extension must be [{'|'.join(format_list)}]"
        parser.error(err_msg)

    args.gbk_resources = gbk_files2gbk_objects(args.gbk_resources)

    return args


if __name__ == "__main__":
    main()
