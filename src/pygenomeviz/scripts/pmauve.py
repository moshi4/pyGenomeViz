import argparse
import csv
import os
import re
import shutil
import subprocess as sp
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

from pygenomeviz import GenomeViz, __version__
from pygenomeviz.scripts import get_argparser, print_args
from pygenomeviz.utils import ColorCycler


def main():
    """Visualization workflow using progressiveMauve"""
    # Get arguments
    args = get_args()
    print_args(args)

    # General options
    seq_files: List[Path] = args.seq_files
    outdir: Path = args.outdir
    format: str = args.format
    reuse: bool = args.reuse
    # Figure appearence options
    fig_width: float = args.fig_width
    fig_track_height: float = args.fig_track_height
    feature_track_ratio: float = args.feature_track_ratio
    link_track_ratio: float = args.link_track_ratio
    tick_track_ratio: float = args.tick_track_ratio
    track_labelsize: int = args.track_labelsize
    tick_labelsize: int = args.tick_labelsize
    normal_link_color: str = args.normal_link_color
    inverted_link_color: str = args.inverted_link_color
    align_type: str = args.align_type
    tick_style: Optional[str] = args.tick_style
    plotstyle: str = args.plotstyle
    cmap: str = args.cmap
    curve: bool = args.curve

    # Setup output contents
    seq_outdir = outdir / "seqfiles"
    os.makedirs(seq_outdir, exist_ok=True)
    result_fig_file = outdir / f"result.{format}"
    xmfa_file = outdir / "mauve.xmfa"
    bbone_file = outdir / "mauve_bbone.tsv"
    ColorCycler.set_cmap(cmap)

    # Copy seqfiles to output directory
    seq_copy_files = []
    for seq_file in seq_files:
        seq_copy_file = seq_outdir / seq_file.name
        shutil.copy(seq_file, seq_copy_file)
        seq_copy_files.append(str(seq_copy_file))

    # Run progressiveMauve
    if reuse and xmfa_file.exists() and bbone_file.exists():
        print("Reuse previous progressiveMauve result.")
    else:
        cmd = f"progressiveMauve --output {xmfa_file} --backbone-output={bbone_file} "
        cmd += f"{' '.join(seq_copy_files)}"
        sp.run(cmd, shell=True)

    # Parse progressiveMauve results
    filenames = parse_xmfa_file(xmfa_file)
    seqid2maxsize, seqid2links = parse_bbone_file(bbone_file)

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

    name2maxsize: Dict[str, int] = defaultdict(int)
    name2blocks: Dict[str, List[Tuple[int, int, int]]] = defaultdict(list)
    for i in range(len(filenames)):
        name2maxsize[filenames[i]] = seqid2maxsize[i]
        name2blocks[filenames[i]] = seqid2links[i]

    # Set tracks & features(blocks)
    for name in filenames:
        maxsize = name2maxsize[name]
        track = gv.add_feature_track(name, maxsize, track_labelsize)
        blocks = name2blocks[name]
        colors = ColorCycler.get_color_list(len(blocks))
        for block, color in zip(blocks, colors):
            start, end, strand = block
            track.add_feature(start, end, strand, plotstyle=plotstyle, facecolor=color)

    # Set links
    for i in range(len(filenames) - 1):
        name1, name2 = filenames[i], filenames[i + 1]
        blocks1, blocks2 = name2blocks[name1], name2blocks[name2]
        for block1, block2 in zip(blocks1, blocks2):
            start1, end1, strand1 = block1
            if strand1 == -1:
                start1, end1 = end1, start1
            start2, end2, strand2 = block2
            if strand2 == -1:
                start2, end2 = end2, start2
            link1, link2 = (name1, start1, end1), (name2, start2, end2)
            gv.add_link(
                link1, link2, normal_link_color, inverted_link_color, curve=curve
            )

    # Save figure
    gv.savefig(result_fig_file, dpi=300)
    print(f"\nSave result figure ({result_fig_file}).")


def parse_xmfa_file(xmfa_file: Union[str, Path]) -> List[str]:
    """Parse progressiveMauve xmfa file

    Parameters
    ----------
    xmfa_file : Union[str, Path]
        progressiveMauve xmfa format file

    Returns
    -------
    filenames : List[str]
        File names
    """
    filenames = []
    with open(xmfa_file) as f:
        lines = f.read().splitlines()
        for line in lines:
            if not line.startswith("#"):
                break
            if re.match("^#Sequence[0-9]+File", line):
                seqfile = Path(line.split("\t")[1])
                filename = seqfile.with_suffix("").name
                filenames.append(filename)
    return filenames


def parse_bbone_file(
    bbone_file: Union[str, Path],
) -> Tuple[Dict[int, int], Dict[int, List[Tuple[int, int, int]]]]:
    """Parse progressiveMauve bbone file

    Parameters
    ----------
    bbone_file : Union[str, Path]
        progressiveMauve bbone format file

    Returns
    -------
    Tuple[Dict[int, int], Dict[int, List[Tuple[int, int, int]]]]
        seqid2maxsize, seqid2blocks
    """
    refid = 0
    with open(bbone_file) as f:
        reader = csv.reader(f, delimiter="\t")
        header_row = next(reader)
        genome_num = int(len(header_row) / 2)
        rows = []
        for row in reader:
            row = [int(col) for col in row]
            # Always set reference seq coordinates to positive value
            if row[refid] < 0:
                row = [col * -1 for col in row]
            # Ignore no commonly conserved regions in all genomes
            if row.count(0) >= 2:
                continue
            rows.append(row)
        # Sort by reference seq coordinates
        rows = sorted(rows, key=lambda row: row[refid])

    seqid2maxsize: Dict[int, int] = defaultdict(int)
    seqid2blocks: Dict[int, List[Tuple[int, int, int]]] = defaultdict(list)
    for row in rows:
        for i in range(genome_num):
            idx = i * 2
            left, right = row[idx], row[idx + 1]
            strand = -1 if left < 0 or right < 0 else 1
            start, end = abs(left), abs(right)
            if seqid2maxsize[i] < end:
                seqid2maxsize[i] = end
            seqid2blocks[i].append((start, end, strand))
    return seqid2maxsize, seqid2blocks


def get_args() -> argparse.Namespace:
    """Get arguments

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
    default_format = "png"
    format_list = ["png", "jpg", "svg", "pdf"]
    general_opts.add_argument(
        "--format",
        type=str,
        help="Output image format ('png'[*]|'jpg'|'svg'|'pdf')",
        default=default_format,
        choices=format_list,
        metavar="",
    )
    general_opts.add_argument(
        "--reuse",
        help="Reuse previous alignment result if available",
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
    args = parser.parse_args()

    # Check arguments
    err_info = ""
    if len(args.seq_files) < 2:
        err_info += "--seq_files must be set at least two files.\n"
    for file in args.seq_files:
        if not file.exists():
            err_info += f"File not found '{file}'.\n"

    if err_info != "":
        parser.error("\n" + err_info)

    return args


if __name__ == "__main__":
    main()
