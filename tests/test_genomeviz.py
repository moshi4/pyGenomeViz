import random
from pathlib import Path

import pytest
from pygenomeviz import Genbank, GenomeViz, load_dataset

random.seed(0)


def test_manual_dataset(tmp_path: Path):
    """Test with simple manual dataset"""
    genome_list = (
        {
            "name": "genome 01",
            "size": 1000,
            "cds_list": ((150, 300, 1), (500, 700, -1), (750, 950, 1)),
        },
        {
            "name": "genome 02",
            "size": 1300,
            "cds_list": ((50, 200, 1), (350, 450, 1), (700, 900, -1), (950, 1150, -1)),
        },
        {
            "name": "genome 03",
            "size": 1200,
            "cds_list": ((150, 300, 1), (350, 450, -1), (500, 700, -1), (701, 900, -1)),
        },
    )

    gv = GenomeViz()
    for genome in genome_list:
        name, size, cds_list = genome["name"], genome["size"], genome["cds_list"]
        track = gv.add_feature_track(name, size)
        for idx, cds in enumerate(cds_list, 1):
            plotstyle = random.choice(("bigarrow", "arrow", "bigbox", "box"))
            labelvpos = random.choice(("top", "center", "bottom", "strand"))
            labelhpos = random.choice(("left", "center", "right"))
            labelha = random.choice(("left", "center", "right"))
            track.add_feature(
                *cds,
                label=f"gene{idx:02d}",
                plotstyle=plotstyle,
                labelvpos=labelvpos,
                labelhpos=labelhpos,
                labelha=labelha,
            )

    gv.add_link(("genome 01", 150, 300), ("genome 02", 50, 200))
    gv.add_link(("genome 01", 700, 500), ("genome 02", 900, 700))
    gv.add_link(("genome 01", 750, 950), ("genome 02", 1150, 950))
    params = {"normal_color": "skyblue", "inverted_color": "lime", "curve": True}
    gv.add_link(("genome 02", 50, 200), ("genome 03", 150, 300), **params)
    gv.add_link(("genome 02", 350, 450), ("genome 03", 450, 350), **params)
    gv.add_link(("genome 02", 900, 700), ("genome 03", 700, 500), **params)
    gv.add_link(("genome 03", 900, 701), ("genome 02", 1150, 950), **params)

    fig_outfile = tmp_path / "manual_dataset.jpg"
    gv.savefig(fig_outfile, dpi=300)
    gv.print_tracks_info(detail=True)

    assert fig_outfile.exists()


def test_escherichia_phage_dataset(tmp_path: Path):
    """Test with 'escherichia phage' dataset"""
    gv = GenomeViz(
        fig_width=12,
        feature_track_ratio=0.5,
        align_type="center",
    )

    gbk_files, links = load_dataset("escherichia_phage")
    for gbk_file in gbk_files:
        gbk = Genbank(gbk_file)
        track = gv.add_feature_track(gbk.name, gbk.genome_length)
        track.add_genbank_features(
            gbk, size_ratio=0.5, arrow_shaft_ratio=1.0, linewidth=0.5
        )

    for link in links:
        link_data1 = (link.ref_name, link.ref_start, link.ref_end)
        link_data2 = (link.query_name, link.query_start, link.query_end)
        gv.add_link(link_data1, link_data2, curve=True)

    fig_outfile = tmp_path / "escherichia_phage.png"
    gv.savefig(fig_outfile, dpi=300)

    assert fig_outfile.exists()


def test_erwinia_phage_dataset(tmp_path: Path):
    """Test with 'erwinia phage' dataset"""
    gv = GenomeViz(
        fig_track_height=0.8,
        feature_track_ratio=0.5,
        tick_track_ratio=0.3,
        tick_style="axis",
    )

    gbk_files, links = load_dataset("erwinia_phage")
    for gbk_file in gbk_files:
        gbk = Genbank(gbk_file)
        track = gv.add_feature_track(gbk.name, gbk.genome_length)
        track.add_genbank_features(gbk, plotstyle="arrow")

    min_identity = int(min(link.identity for link in links))
    for link in links:
        link_data1 = (link.ref_name, link.ref_start, link.ref_end)
        link_data2 = (link.query_name, link.query_start, link.query_end)
        gv.add_link(link_data1, link_data2, v=link.identity, vmin=min_identity)

    fig = gv.plotfig()
    gv.set_colorbar(fig, vmin=min_identity)
    fig_outfile = tmp_path / "erwinia_phage.pdf"
    fig.savefig(fig_outfile)

    assert fig_outfile.exists()


def test_enterobacteria_phage_dataset(tmp_path: Path):
    """Test with 'enterobacteria phage' dataset"""
    gv = GenomeViz(
        fig_width=10,
        fig_track_height=0.7,
        feature_track_ratio=0.5,
        tick_track_ratio=0.5,
        align_type="center",
        tick_style="bar",
        tick_labelsize=10,
    )

    gbk_files, links = load_dataset("enterobacteria_phage")
    for idx, gbk_file in enumerate(gbk_files):
        gbk = Genbank(gbk_file)
        track = gv.add_feature_track(gbk.name, gbk.genome_length, labelsize=10)
        track.add_genbank_features(
            gbk,
            label_type="product" if idx == 0 else None,  # Labeling only top track
            label_filter=["hypothetical"],  # Ignore 'hypothetical ~~~' label
            labelsize=8,
            labelvpos="top",
            facecolor="skyblue",
            linewidth=0.5,
        )

    normal_color, inverted_color, alpha = "chocolate", "limegreen", 0.5
    min_identity = int(min(link.identity for link in links))
    for link in links:
        link_data1 = (link.ref_name, link.ref_start, link.ref_end)
        link_data2 = (link.query_name, link.query_start, link.query_end)
        gv.add_link(
            link_data1,
            link_data2,
            normal_color,
            inverted_color,
            alpha,
            v=link.identity,
            vmin=min_identity,
            curve=True,
        )

    fig = gv.plotfig()
    gv.set_colorbar(
        fig,
        bar_colors=[normal_color, inverted_color],
        alpha=alpha,
        vmin=min_identity,
        bar_height=0.15,
        bar_label="Identity",
        bar_labelsize=10,
    )
    fig_outfile = tmp_path / "enterobacteria_phage.svg"
    fig.savefig(fig_outfile, dpi=300, bbox_inches="tight", pad_inches=0.5)


def test_top_track():
    """Test top track"""
    gv = GenomeViz()
    top_track_name, top_track_size = "Top Track", 1000
    gv.add_feature_track(top_track_name, top_track_size)
    gv.add_feature_track("track02", 1200)

    assert gv.top_track.name == top_track_name
    assert gv.top_track.size == top_track_size


def test_add_feature_subtrack():
    """Test add feature subtrack"""
    gv = GenomeViz()
    top_track_name, top_track_size = "Top Track", 1000
    gv.add_feature_track(top_track_name, top_track_size)
    gv.add_feature_track("track02", 1200)

    gv.add_feature_subtrack(top_track_name, "subtrack01")
    gv.add_feature_subtrack(top_track_name, "subtrack02")

    assert len(gv.top_track.subtracks) == 2
    assert gv.top_track.subtracks[0].name == "subtrack01"
    assert gv.top_track.subtracks[1].name == "subtrack02"


def test_ax_property():
    """Test ax property access"""
    gv = GenomeViz()
    gv.add_feature_track("track01", 1000)
    # Can't access ax property before calling 'plotfig' method
    with pytest.raises(ValueError):
        gv.top_track.ax
    gv.plotfig()
    gv.top_track.ax


def test_offset_property():
    """Test offset property access"""
    gv = GenomeViz()
    gv.add_feature_track("track01", 1000)
    # Can't access offset property before calling 'plotfig' method
    with pytest.raises(ValueError):
        gv.top_track.offset
    gv.plotfig()
    gv.top_track.offset
