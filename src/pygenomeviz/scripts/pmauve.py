#!/usr/bin/env python
from __future__ import annotations

import argparse
import os
import signal
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Sequence

from pygenomeviz import GenomeViz
from pygenomeviz.align import AlignCoord, ProgressiveMauve
from pygenomeviz.logger import get_logger
from pygenomeviz.scripts import (
    ALIGN_COORDS_FILENAME,
    LOG_FILENAME,
    CustomHelpFormatter,
    log_basic_env_info,
    setup_argparser,
    validate_args,
)
from pygenomeviz.typing import PlotStyle, TrackAlignType
from pygenomeviz.utils import ColorCycler

CLI_NAME = "pgv-pmauve"


def main():
    """Main function called from CLI"""
    args = get_args()
    args_dict = args.__dict__
    try:
        run(**args_dict, log_params=args_dict)
    except KeyboardInterrupt:
        get_logger(__name__).error("Keyboard Interrupt")
        sys.exit(-signal.SIGINT)
    except Exception as e:
        get_logger(__name__).error(e)
        sys.exit(getattr(e, "errno", 1))


def run(
    # General options
    seqs: Sequence[str | Path],
    outdir: str | Path,
    formats: list[str],
    reuse: bool,
    quiet: bool,
    debug: bool,
    # Figure appearence options
    fig_width: float,
    fig_track_height: float,
    track_align_type: TrackAlignType,
    feature_track_ratio: float,
    show_scale_bar: bool,
    show_scale_xticks: bool,
    curve: bool,
    dpi: int,
    track_labelsize: int,
    scale_labelsize: int,
    normal_link_color: str,
    inverted_link_color: str,
    refid: int,
    block_plotstyle: PlotStyle,
    block_cmap: str,
    # Log parameters
    log_params: dict | None = None,
):
    """Run genome visualization workflow"""
    start_time = time.time()

    # Make output directory
    outdir = Path(outdir)
    os.makedirs(outdir, exist_ok=True)

    # Set logger
    log_file = outdir / LOG_FILENAME
    logger = get_logger(__name__, log_file, quiet)
    log_basic_env_info(logger, CLI_NAME, log_params)

    # Run progressiveMauve alignment
    align_coords_file = outdir / ALIGN_COORDS_FILENAME
    pmauve = ProgressiveMauve(
        seqs,
        outdir=outdir / "tmp" if debug else None,
        refid=refid,
        logger=logger,
        quiet=quiet,
    )
    if reuse and align_coords_file.exists():
        logger.info(f"Reuse alignment result in '{align_coords_file}'")
        align_coords = AlignCoord.read(align_coords_file)
    else:
        align_coords = pmauve.run()
        logger.info(f"Write alignment result to '{align_coords_file}'")
        AlignCoord.write(align_coords, align_coords_file)

    # Setup synteny blocks
    name2blocks: dict[str, list[tuple[int, int, int]]] = defaultdict(list)
    for ac in align_coords:
        if ac.query_block not in name2blocks[ac.query_name]:
            name2blocks[ac.query_name].append(ac.query_block)
        if ac.ref_block not in name2blocks[ac.ref_name]:
            name2blocks[ac.ref_name].append(ac.ref_block)

    # Create GenomeViz instance
    gv = GenomeViz(
        fig_width=fig_width,
        fig_track_height=fig_track_height,
        track_align_type=track_align_type,
        feature_track_ratio=feature_track_ratio,
    )
    if show_scale_bar:
        gv.set_scale_bar(labelsize=scale_labelsize)
    if show_scale_xticks:
        gv.set_scale_xticks(labelsize=scale_labelsize)

    ColorCycler.set_cmap(block_cmap)
    for name, seqlen in pmauve.name2seqlen.items():
        # Add track
        track = gv.add_feature_track(
            name,
            {name: seqlen},
            labelsize=track_labelsize,
        )
        # Add synteny blocks
        blocks = name2blocks[name]
        colors = ColorCycler.get_color_list(len(blocks))
        for block, color in zip(blocks, colors):
            track.add_feature(*block, plotstyle=block_plotstyle, fc=color)

    # Plot links
    for ac in align_coords:
        gv.add_link(
            ac.query_link,
            ac.ref_link,
            curve=curve,
            color=normal_link_color,
            inverted_color=inverted_link_color,
        )

    # Output result image file
    for format in formats:
        output_file = outdir / f"result.{format}"
        logger.info(f"Plotting {format.upper()} format result...")
        if format == "html":
            gv.savefig_html(output_file)
        else:
            gv.savefig(output_file, dpi=dpi)
        logger.info(f"Output result image file '{output_file}'")

    elapsed_time = time.time() - start_time
    logger.info(f"Done (elapsed time: {elapsed_time:.2f}[s])")


def get_args() -> argparse.Namespace:
    """Get arguments

    Returns
    -------
    args : argparse.Namespace
        Argument parameters
    """
    parser = argparse.ArgumentParser(
        description="pyGenomeViz CLI workflow using progressiveMauve",
        usage=f"{CLI_NAME} [options] seq1.gbk seq2.gbk seq3.gbk -o outdir",
        epilog="[*] marker means the default value.",
        add_help=False,
        allow_abbrev=False,
        formatter_class=CustomHelpFormatter,
    )
    setup_argparser(parser, CLI_NAME)

    args = parser.parse_args()
    validate_args(args, parser)

    return args


if __name__ == "__main__":
    main()
