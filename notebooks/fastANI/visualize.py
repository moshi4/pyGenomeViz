#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

from Bio import SeqIO
from pygenomeviz import GenomeViz
from pygenomeviz.utils import ColorCycler


def main():
    """Run visualize conserved regions workflow"""
    # Parse arguments
    args = get_args()
    genome_fasta_file1: Path = args.fasta_file1
    genome_fasta_file2: Path = args.fasta_file2
    fastani_visual_file: Path = args.visual_file
    visualize_result_file: Path = args.plot_outfile
    cmap: str = args.cmap
    link_color: str = args.link_color
    curve: bool = args.curve

    # Load genome fasta information
    genome_name1 = genome_fasta_file1.with_suffix("").name
    genome_name2 = genome_fasta_file2.with_suffix("").name
    records1 = SeqIO.parse(genome_fasta_file1, "fasta")
    seq_length1 = sum([len(r) for r in records1])
    records2 = SeqIO.parse(genome_fasta_file2, "fasta")
    seq_length2 = sum([len(r) for r in records2])

    # Load fastANI visual result
    fastani_results = []
    with open(fastani_visual_file) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            # Ignore blank lines
            if len(row) == 0:
                continue
            start1, end1 = int(row[6]), int(row[7])
            start2, end2 = int(row[8]), int(row[9])
            identity = float(row[2])
            link1, link2 = (genome_name1, start1, end1), (genome_name2, start2, end2)
            fastani_results.append((link1, link2, identity))

    # Visualize conserved regions detected by fastANI
    gv = GenomeViz(
        fig_width=15,
        fig_track_height=1.0,
        feature_track_ratio=0.1,
        tick_track_ratio=0.2,
        align_type="center",  # "left", "center", "right"
        tick_style="bar",  # "axis", "bar", None
        plot_size_thr=0,
    )
    track1 = gv.add_feature_track(genome_name1, seq_length1)
    track2 = gv.add_feature_track(genome_name2, seq_length2)

    ColorCycler.set_cmap(cmap)  # "hsv", "viridis", "jet", etc...
    colors = ColorCycler.get_color_list(len(fastani_results))

    min_identity = int(min([res[2] for res in fastani_results]))
    for res, color in zip(fastani_results, colors):
        link1, link2, identity = res
        track1.add_feature(link1[1], link1[2], plotstyle="bigbox", facecolor=color)
        track2.add_feature(link2[1], link2[2], plotstyle="bigbox", facecolor=color)
        gv.add_link(
            link1, link2, link_color, v=identity, vmin=min_identity, curve=curve
        )

    fig = gv.plotfig()
    gv.set_colorbar(
        fig, bar_colors=[link_color], vmin=min_identity, bar_height=0.3, bar_bottom=0.15
    )
    fig.savefig(visualize_result_file)


def get_args() -> argparse.Namespace:
    """Get arguments
    Returns:
        argparse.Namespace: Argument values
    """
    description = "Visualize conserved regions detected by fastANI"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "fasta_file1",
        type=Path,
        help="Input genome fasta 1",
    )
    parser.add_argument(
        "fasta_file2",
        type=Path,
        help="Input genome fasta 2",
    )
    parser.add_argument(
        "visual_file",
        type=Path,
        help="fastANI visual result file",
    )
    parser.add_argument(
        "plot_outfile",
        type=Path,
        help="Plot result outfile [*.jpg|*.png|*.svg|*.pdf]",
    )
    default_cmap = "hsv"
    parser.add_argument(
        "--cmap",
        type=str,
        help=f"Colormap for matplotlib (Default: '{default_cmap}')",
        default=default_cmap,
        metavar="",
    )
    default_link_color = "grey"
    parser.add_argument(
        "--link_color",
        type=str,
        help=f"Link color (Default: '{default_link_color}')",
        default=default_link_color,
        metavar="",
    )
    parser.add_argument(
        "--curve",
        help="Plot curved style link (Default: OFF)",
        action="store_true",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
