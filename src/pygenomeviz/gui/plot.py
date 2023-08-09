from __future__ import annotations

import hashlib
import os
from pathlib import Path
from tempfile import TemporaryDirectory

from matplotlib.figure import Figure

from pygenomeviz import Genbank, GenomeViz
from pygenomeviz.align import AlignCoord, MMseqs, MUMmer
from pygenomeviz.gui import config, utils


def create_genomeviz(
    gbk_list: list[Genbank],
    cfg: config.PgvConfig,
) -> tuple[GenomeViz, Figure, list[AlignCoord]]:
    """Create GenomeViz from genbank list

    Parameters
    ----------
    gbk_list : list[Genbank]
        Genbank list
    cfg : PgvConfig
        Config

    Returns
    -------
    gv, fig : tuple[GenomeViz, Figure]
        GenomeViz, Figure
    """
    # Create genomeviz instance
    gv = GenomeViz(
        fig_width=cfg.fig.width,
        fig_track_height=cfg.fig.track_height,
        feature_track_ratio=cfg.fig.track_ratio,
        link_track_ratio=1.0,
        tick_track_ratio=0.5,
        align_type=cfg.fig.align_type,  # type: ignore
        tick_style=cfg.fig.tick_style,  # type: ignore
    )

    # Add features from genbank file
    for gbk_cnt, gbk in enumerate(gbk_list):
        track = gv.add_feature_track(
            name=gbk.name,
            size=gbk.range_size,
            start_pos=gbk.min_range,
            labelsize=cfg.fig.label_size,
        )
        track.set_sublabel(
            size=cfg.fig.range_label_size,
            sublabel_kws=dict(
                bbox=dict(fc="white", ec="none", alpha=0.5, boxstyle="square,pad=0.0")
            ),
        )
        if cfg.feat.show_only_top_label and gbk_cnt != 0:
            label_type = None
        else:
            label_type = cfg.feat.label_type

        for type in cfg.feat.types:
            if type == "Pseudo":
                feature_type, pseudogene = "CDS", True
            else:
                feature_type, pseudogene = type, False
            track.add_genbank_features(
                gbk,
                feature_type=feature_type,
                pseudogene=pseudogene,
                plotstyle=cfg.feat.type2plotstyle[type],  # type: ignore
                facecolor=cfg.feat.type2color[type],
                label_handle_func=cfg.feat.label_filter_func,
                linewidth=0.5,
                arrow_shaft_ratio=0.5,
                labelsize=cfg.feat.label_size,
                label_type=label_type,
                labelvpos="top",
            )

    # Return if alignment process is not required
    if cfg.aln.method is None or len(gbk_list) == 1:
        return gv, gv.plotfig(), []

    # Create processig cache directory
    package_name = __name__.split(".")[0]
    gui_cache_dir = Path.home() / ".cache" / package_name / "gui"
    os.makedirs(gui_cache_dir, exist_ok=True)
    utils.remove_old_files(gui_cache_dir)

    # Create md5 hash unique filename to enable alignment result cache
    md5_hash_source = "\n".join([str(gbk) for gbk in gbk_list]).encode()
    md5_hash_value = hashlib.md5(md5_hash_source).hexdigest()
    aln_coords_filename = f"{md5_hash_value}_{cfg.aln.method}.tsv"
    aln_coords_file = gui_cache_dir / aln_coords_filename.replace(" ", "")

    # Genome alignment
    if aln_coords_file.exists():
        align_coords = AlignCoord.read(aln_coords_file)
    else:
        with TemporaryDirectory() as tmpdir:
            if cfg.aln.method == "MUMmer (protein)":
                aligner = MUMmer(gbk_list, tmpdir, "protein", "many-to-many", 1)
            elif cfg.aln.method == "MUMmer (nucleotide)":
                aligner = MUMmer(gbk_list, tmpdir, "nucleotide", "many-to-many", 1)
            elif cfg.aln.method == "MMseqs":
                aligner = MMseqs(gbk_list, tmpdir, quiet=True)
            else:
                raise ValueError(f"{cfg.aln.method=} is invalid.")
            align_coords = aligner.run()
            AlignCoord.write(align_coords, aln_coords_file)
    align_coords = AlignCoord.filter(
        align_coords, cfg.aln.min_length, cfg.aln.min_identity
    )

    if len(align_coords) == 0:
        return gv, gv.plotfig(), []

    # Add alignment links
    min_identity = int(min([ac.identity for ac in align_coords]))
    for ac in align_coords:
        gv.add_link(
            ac.ref_link,
            ac.query_link,
            cfg.aln.normal_link_color,
            cfg.aln.inverted_link_color,
            curve=cfg.aln.curve,
            v=ac.identity,
            vmin=min_identity,
        )

    # Plot figure
    fig = gv.plotfig()

    # Set colorbar
    bar_colors = [cfg.aln.normal_link_color]
    has_inverted_link = any([ac.is_inverted for ac in align_coords])
    if has_inverted_link:
        bar_colors.append(cfg.aln.inverted_link_color)
    gv.set_colorbar(
        fig,
        bar_colors,
        vmin=min_identity,
        bar_height=cfg.aln.colorbar_height,
        tick_labelsize=15,
    )

    return gv, fig, align_coords
