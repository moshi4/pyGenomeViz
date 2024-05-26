from __future__ import annotations

import hashlib
import os
from pathlib import Path

from pygenomeviz import GenomeViz
from pygenomeviz.align import AlignCoord, AlignToolBase, Blast, MMseqs, MUMmer
from pygenomeviz.exception import SegmentNotFoundError
from pygenomeviz.gui import config, utils
from pygenomeviz.parser import Genbank
from pygenomeviz.typing import AlnMethod
from pygenomeviz.utils import is_pseudo_feature


def plot_by_gui_cfg(
    gbk_list: list[Genbank],
    cfg: config.PgvGuiPlotConfig,
) -> tuple[GenomeViz, list[AlignCoord]]:
    """Plot by GUI configs

    Parameters
    ----------
    gbk_list : list[Genbank]
        Genbank list
    cfg : PgvConfig
        Config

    Returns
    -------
    gv : GenomeViz
        GenomeViz instance
    align_coords : list[AlignCoord]
        AlignCoord list
    """
    # Create genomeviz instance
    gv = GenomeViz(
        fig_width=cfg.fig.width,
        fig_track_height=cfg.fig.track_height,
        track_align_type=cfg.fig.track_align_type,  # type: ignore
        feature_track_ratio=cfg.fig.feature_track_ratio,
        link_track_ratio=cfg.fig.link_track_ratio,
    )
    if cfg.fig.scale_style == "bar":
        gv.set_scale_bar()
    elif cfg.fig.scale_style == "xticks":
        gv.set_scale_xticks()

    # Add features from genbank file
    for gbk_cnt, gbk in enumerate(gbk_list):
        seqid2range = cfg.name2seqid2range[gbk.name]
        track = gv.add_feature_track(
            name=gbk.name,
            segments=seqid2range,
            space=cfg.fig.seg_space_ratio,
            labelsize=cfg.fig.label_size,
        )
        for seg in track.segments:
            seg.add_sublabel(
                size=cfg.fig.range_label_size,
                bbox=dict(fc="white", ec="none", alpha=0.5, boxstyle="square,pad=0.0"),
            )
        label_type = cfg.feat.label_type
        if cfg.feat.label_target_track == "top" and gbk_cnt != 0:
            label_type = None
        seqid2features = gbk.get_seqid2features(feature_type=cfg.feat.types)
        for seqid, _ in seqid2range.items():
            features = seqid2features[seqid]
            for feature in features:
                fc = cfg.feat.type2color[feature.type]
                if is_pseudo_feature(feature):
                    fc = cfg.feat.pseudo_color
                track.add_features(
                    feature,
                    target_seg=seqid,
                    plotstyle=cfg.feat.type2plotstyle[feature.type],  # type: ignore
                    arrow_shaft_ratio=0.5,
                    label_type=label_type,
                    label_handler=cfg.feat.label_filter_func,
                    ignore_outside_range=True,
                    text_kws=dict(size=cfg.feat.label_size, vpos="top"),
                    fc=fc,
                    lw=cfg.feat.line_width,
                )

    # Return if alignment process is not required
    if cfg.aln.method is None or len(gbk_list) == 1:
        return gv, []

    # Create processig cache directory
    package_name = __name__.split(".")[0]
    gui_cache_dir = Path.home() / ".cache" / package_name / "gui"
    os.makedirs(gui_cache_dir, exist_ok=True)
    utils.remove_old_files(gui_cache_dir)

    # Create md5 hash unique filename to enable cache alignment result
    md5_hash_source = "\n".join([f"{gbk} {gbk.full_genome_seq}" for gbk in gbk_list])
    md5_hash_value = hashlib.md5(md5_hash_source.encode()).hexdigest()
    aln_coords_filename = f"{md5_hash_value}_{cfg.aln.method}.tsv"
    aln_coords_file = gui_cache_dir / aln_coords_filename.replace(" ", "")

    # Genome alignment
    if aln_coords_file.exists():
        align_coords = AlignCoord.read(aln_coords_file)
    else:
        aln_method2aligner: dict[AlnMethod, AlignToolBase] = {
            "MUMmer (nucleotide)": MUMmer(gbk_list, seqtype="nucleotide", quiet=True),
            "MUMmer (protein)": MUMmer(gbk_list, seqtype="protein", quiet=True),
            "MMseqs RBH": MMseqs(gbk_list, quiet=True),
            "BLAST (nucleotide)": Blast(gbk_list, seqtype="nucleotide", quiet=True),
            "BLAST (protein)": Blast(gbk_list, seqtype="protein", quiet=True),
        }
        aligner = aln_method2aligner[cfg.aln.method]
        align_coords = aligner.run()
        AlignCoord.write(align_coords, aln_coords_file)
    align_coords = AlignCoord.filter(
        align_coords,
        length_thr=cfg.aln.min_length,
        identity_thr=cfg.aln.min_identity,
    )

    if len(align_coords) == 0:
        return gv, []

    # Add alignment links
    min_identity = int(min([ac.identity for ac in align_coords if ac.identity]))
    for ac in align_coords:
        try:
            gv.add_link(
                (ac.query_id, ac.query_name, ac.query_start, ac.query_end),
                (ac.ref_id, ac.ref_name, ac.ref_start, ac.ref_end),
                color=cfg.aln.normal_link_color,
                inverted_color=cfg.aln.inverted_link_color,
                curve=cfg.aln.curve,
                v=ac.identity,
                vmin=min_identity,
                ignore_outside_range=True,
            )
        except SegmentNotFoundError:
            continue

    # Set colorbar
    bar_colors = [cfg.aln.normal_link_color]
    has_inverted_link = any([ac.is_inverted for ac in align_coords])
    if has_inverted_link:
        bar_colors.append(cfg.aln.inverted_link_color)
    gv.set_colorbar(
        bar_colors,
        vmin=min_identity,
        bar_height=cfg.aln.colorbar_height,
    )

    return gv, align_coords
