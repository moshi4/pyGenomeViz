import argparse
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import List, Optional

from pygenomeviz import Genbank, GenomeViz, __version__
from pygenomeviz.align import Align, AlignCoord


def main():
    """Visualization workflow using progressiveMauve"""
    # Get arguments
    args = get_args()
    gbk_files: List[Path] = args.gbk_files
    outdir: Path = args.outdir
    format: str = args.format
    fig_width: float = args.fig_width
    fig_track_height: float = args.fig_track_height
    feature_track_ratio: float = args.feature_track_ratio
    link_track_ratio: float = args.link_track_ratio
    tick_track_ratio: float = args.tick_track_ratio
    track_labelsize: int = args.track_labelsize
    tick_labelsize: int = args.tick_labelsize
    normal_link_color: str = args.normal_link_color
    inverted_link_color: str = args.inverted_link_color
    curve: bool = args.curve
    align_type: str = args.align_type
    tick_style: Optional[str] = args.tick_style
    plotstyle: str = args.plotstyle
    seqtype: str = args.seqtype
    maptype: str = args.maptype
    force: bool = args.force

    os.makedirs(outdir, exist_ok=True)
    result_file = outdir / f"result.{format}"
    align_coords_file = outdir / "align_coords.tsv"

    gv = GenomeViz(
        fig_width=fig_width,
        fig_track_height=fig_track_height,
        feature_track_ratio=feature_track_ratio,
        link_track_ratio=link_track_ratio,
        tick_track_ratio=tick_track_ratio,
        align_type=align_type,
        tick_style=tick_style,
        tick_labelsize=tick_labelsize,
        plot_size_thr=0.0005,
    )

    gbk_list = [Genbank(f) for f in gbk_files]
    for gbk in gbk_list:
        track = gv.add_feature_track(gbk.name, gbk.genome_length, track_labelsize)
        track.add_genbank_features(
            gbk, plotstyle=plotstyle, facecolor="orange", linewidth=0.0
        )

    if force or not align_coords_file.exists():
        # Run MUMmer alignment
        with TemporaryDirectory() as tmpdir:
            align_coords = Align(gbk_list, tmpdir, seqtype, maptype).run()
            AlignCoord.write(align_coords, align_coords_file)
    else:
        align_coords = AlignCoord.read(align_coords_file)
    # Filter alignment results by 'min_length', 'min_identity' threshold
    align_coords = AlignCoord.filter(align_coords, min_length=0, min_identity=0)

    min_identity = int(min([ac.identity for ac in align_coords]))
    for ac in align_coords:
        link1 = (ac.ref_name, ac.ref_start, ac.ref_end)
        link2 = (ac.query_name, ac.query_start, ac.query_end)
        gv.add_link(
            link1,
            link2,
            normal_link_color,
            inverted_link_color,
            curve=curve,
            v=ac.identity,
            vmin=min_identity,
        )

    gv.savefig(result_file)


def get_args() -> argparse.Namespace:
    """Get arguments

    Returns
    -------
    args : argparse.Namespace
        Argument parameters
    """

    class CustomHelpFormatter(argparse.HelpFormatter):
        def __init__(self, prog, indent_increment=2, max_help_position=30, width=None):
            super().__init__(prog, indent_increment, max_help_position, width)

    description = "pyGenomeViz visualization workflow using progressiveMauve"
    parser = argparse.ArgumentParser(
        description=description, add_help=False, formatter_class=CustomHelpFormatter
    )

    parser.add_argument(
        "--gbk_files",
        type=Path,
        help="Input genome genbank files",
        nargs="+",
        required=True,
        metavar="IN",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        help="Output directory",
        required=True,
        metavar="OUT",
    )
    default_format = "png"
    format_list = ["png", "jpg", "svg", "pdf"]
    parser.add_argument(
        "--format",
        type=str,
        help="Output image format ('png'[default]|'jpg'|'svg'|'pdf')",
        default=default_format,
        choices=format_list,
        metavar="",
    )
    defaullt_seqtype = "nucleotide"
    parser.add_argument(
        "--seqtype",
        type=str,
        help="MUMmer alignment sequence type ('nucleotide'[default]|'protein')",
        default=defaullt_seqtype,
        choices=["nucleotide", "protein"],
        metavar="",
    )
    default_maptype = "one-to-one"
    parser.add_argument(
        "--maptype",
        type=str,
        help="MUMmer alignment map type ('one-to-one'[default]|'many-to-many')",
        default=default_maptype,
        choices=["one-to-one", "many-to-many"],
        metavar="",
    )
    default_fig_width = 15
    parser.add_argument(
        "--fig_width",
        type=float,
        help=f"Figure width (Default: {default_fig_width})",
        default=default_fig_width,
        metavar="",
    )
    default_fig_track_height = 1.0
    parser.add_argument(
        "--fig_track_height",
        type=float,
        help=f"Figure track height (Default: {default_fig_track_height})",
        default=default_fig_track_height,
        metavar="",
    )
    default_feature_track_ratio = 1.0
    parser.add_argument(
        "--feature_track_ratio",
        type=float,
        help=f"Feature track ratio (Default: {default_feature_track_ratio})",
        default=default_feature_track_ratio,
        metavar="",
    )
    default_link_track_ratio = 5.0
    parser.add_argument(
        "--link_track_ratio",
        type=float,
        help=f"Link track ratio (Default: {default_link_track_ratio})",
        default=default_link_track_ratio,
        metavar="",
    )
    default_tick_track_ratio = 1.0
    parser.add_argument(
        "--tick_track_ratio",
        type=float,
        help=f"Tick track ratio (Default: {default_tick_track_ratio})",
        default=default_tick_track_ratio,
        metavar="",
    )
    default_track_labelsize = 20
    parser.add_argument(
        "--track_labelsize",
        type=int,
        help=f"Track label size (Default: {default_track_labelsize})",
        default=default_track_labelsize,
        metavar="",
    )
    default_tick_labelsize = 15
    parser.add_argument(
        "--tick_labelsize",
        type=int,
        help=f"Tick label size (Default: {default_tick_labelsize})",
        default=default_tick_labelsize,
        metavar="",
    )
    default_normal_link_color = "grey"
    parser.add_argument(
        "--normal_link_color",
        type=str,
        help=f"Normal link color (Default: '{default_normal_link_color}')",
        default=default_normal_link_color,
        metavar="",
    )
    default_inverted_link_color = "tomato"
    parser.add_argument(
        "--inverted_link_color",
        type=str,
        help=f"Inverted link color (Default: '{default_inverted_link_color}')",
        default=default_inverted_link_color,
        metavar="",
    )
    default_align_type = "center"
    align_type_list = ["left", "center", "right"]
    parser.add_argument(
        "--align_type",
        type=str,
        help="Figure tracks align type ('left'|'center'|'right')",
        default=default_align_type,
        choices=align_type_list,
        metavar="",
    )
    default_tick_style = None
    parser.add_argument(
        "--tick_style",
        type=str,
        help="Tick style ('bar'|'axis'|None[default])",
        default=default_tick_style,
        choices=["bar", "axis"],
        metavar="",
    )
    plotstyle = "bigarrow"
    plotstyle_list = ["bigarrow", "arrow"]
    parser.add_argument(
        "--plotstyle",
        type=str,
        help="Plot feature style ('bigarrow'[default]|'arrow')",
        default=plotstyle,
        choices=plotstyle_list,
        metavar="",
    )
    parser.add_argument(
        "--curve",
        help="Plot curved style link (Default: OFF)",
        action="store_true",
    )
    parser.add_argument(
        "--force",
        help="Forcibly overwrite (not reuse) previous calculation results",
        action="store_true",
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
    args = parser.parse_args()

    # Check arguments
    err_info = ""
    if len(args.gbk_files) < 2:
        err_info += "--seq_files must be set at least two files.\n"
    for file in args.gbk_files:
        if not file.exists():
            err_info += f"File not found '{file}'.\n"

    if err_info != "":
        parser.error("\n" + err_info)

    return args


if __name__ == "__main__":
    main()
