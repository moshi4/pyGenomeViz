#!/usr/bin/env python
from __future__ import annotations

import argparse
import os
import signal
import sys
import time
from pathlib import Path
from typing import Sequence

from pygenomeviz import GenomeViz
from pygenomeviz.align import AlignCoord, Blast
from pygenomeviz.logger import get_logger
from pygenomeviz.parser import Genbank
from pygenomeviz.scripts import (
    ALIGN_COORDS_FILENAME,
    LOG_FILENAME,
    CustomHelpFormatter,
    log_basic_env_info,
    setup_argparser,
    validate_args,
)
from pygenomeviz.typing import PlotStyle, SeqType, TrackAlignType
from pygenomeviz.utils import is_pseudo_feature

CLI_NAME = "pgv-blast"


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
    # Blast alignment options
    seqtype: SeqType,
    threads: int,
    length_thr: int,
    identity_thr: float,
    evalue_thr: float,
    # Figure appearence options
    fig_width: float,
    fig_track_height: float,
    track_align_type: TrackAlignType,
    feature_track_ratio: float,
    show_scale_bar: bool,
    show_scale_xticks: bool,
    track_labelsize: int,
    scale_labelsize: int,
    normal_link_color: str,
    inverted_link_color: str,
    curve: bool,
    dpi: int,
    segment_space: float,
    feature_type2color: dict[str, str],
    pseudo_color: str,
    feature_plotstyle: PlotStyle,
    feature_linewidth: float,
    feature_labeltrack: str,
    feature_labeltype: str | None,
    feature_labelsize: int,
    cbar_width: float,
    cbar_height: float,
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

    # Run Blast alignment
    gbk_list = [Genbank(seq) for seq in seqs]

    align_coords_file = outdir / ALIGN_COORDS_FILENAME
    if reuse and align_coords_file.exists():
        logger.info(f"Reuse alignment result in '{align_coords_file}'")
        align_coords = AlignCoord.read(align_coords_file)
    else:
        align_coords = Blast(
            gbk_list,
            outdir=outdir / "tmp" if debug else None,
            seqtype=seqtype,
            threads=threads,
            logger=logger,
            quiet=quiet,
        ).run()
        logger.info(f"Write alignment result to '{align_coords_file}'")
        AlignCoord.write(align_coords, align_coords_file)
    align_coords = AlignCoord.filter(
        align_coords,
        length_thr=length_thr,
        identity_thr=identity_thr,
        evalue_thr=evalue_thr,
    )
    logger.info(
        f"Filter alignment result by {length_thr=}, {identity_thr=} {evalue_thr=}"
    )

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

    for idx, gbk in enumerate(gbk_list):
        # Add track
        track = gv.add_feature_track(
            gbk.name,
            gbk.get_seqid2size(),
            space=segment_space,
            labelsize=track_labelsize,
        )
        # Add features
        feature_label_type = None
        if feature_labeltrack == "all" or (feature_labeltrack == "top" and idx == 0):
            feature_label_type = feature_labeltype
        seqid2features = gbk.get_seqid2features(feature_type=None)
        for seqid, features in seqid2features.items():
            for feature in features:
                if feature.type in feature_type2color:
                    fc = feature_type2color[feature.type]
                    if is_pseudo_feature(feature):
                        fc = pseudo_color
                    track.add_features(
                        feature,
                        target_seg=seqid,
                        plotstyle=feature_plotstyle,
                        label_type=feature_label_type,
                        fc=fc,
                        lw=feature_linewidth,
                        ec="black",
                        text_kws=dict(size=feature_labelsize),
                    )

    if len(align_coords) > 0:
        min_ident = int(min([ac.identity for ac in align_coords if ac.identity]))
        # Plot links
        for ac in align_coords:
            gv.add_link(
                ac.ref_link,
                ac.query_link,
                curve=curve,
                color=normal_link_color,
                inverted_color=inverted_link_color,
                v=ac.identity,
                vmin=min_ident,
            )
        # Plot colorbar
        colors = [normal_link_color]
        if any([ac.is_inverted for ac in align_coords]):
            colors.append(inverted_link_color)
        gv.set_colorbar(
            colors=colors,
            vmin=min_ident,
            bar_width=cbar_width,
            bar_height=cbar_height,
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
        description="pyGenomeViz CLI workflow using BLAST (blastn, tblastx)",
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
