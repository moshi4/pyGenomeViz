import random
from pathlib import Path

import pytest

from pygenomeviz import Genbank, GenomeViz, Gff, load_dataset

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
    cnt = 0
    for genome in genome_list:
        name, size, cds_list = genome["name"], genome["size"], genome["cds_list"]
        track = gv.add_feature_track(name, size)
        plotstyles = ("bigarrow", "arrow", "bigbox", "box", "bigrbox", "rbox")
        labelvpos_list = ("top", "center", "bottom", "strand")
        labelhpos_list = ("left", "center", "right")
        labelha_list = ("left", "center", "right")
        for cds in cds_list:
            plotstyle = plotstyles[cnt % len(plotstyles)]
            labelvpos = labelvpos_list[cnt % len(labelvpos_list)]
            labelhpos = labelhpos_list[cnt % len(labelhpos_list)]
            labelha = labelha_list[cnt % len(labelha_list)]
            track.add_feature(
                *cds,
                label=f"gene{cnt:02d}",
                plotstyle=plotstyle,
                labelvpos=labelvpos,
                labelhpos=labelhpos,
                labelha=labelha,
            )
            cnt += 1

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
            label_handle_func=lambda s: "" if s.startswith("hypothetical") else s,
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
    assert fig_outfile.exists()


def test_add_exon_feature(tmp_path: Path):
    """Test add exon feature"""
    exon_regions = [(1, 210), (301, 480), (591, 800), (851, 1000), (1031, 1300)]
    exon_labels = [str(i) for i in range(len(exon_regions))]

    gv = GenomeViz(fig_track_height=0.3, link_track_ratio=0.3)
    plotstyles = ("bigarrow", "arrow", "bigbox", "box", "bigrbox", "rbox")
    for plotstyle in plotstyles:
        for strand in (1, -1):
            size = exon_regions[-1][-1]
            track = gv.add_feature_track(name=f"{plotstyle} ({strand})", size=size)
            track.add_exon_feature(
                exon_regions, strand, plotstyle=plotstyle, exon_labels=exon_labels
            )
    fig_outfile = tmp_path / "exon_feature.svg"
    gv.savefig(fig_outfile)
    assert fig_outfile.exists()


def test_gff_genomeviz(gff_file: Path, tmp_path: Path):
    """Test add_gff_features"""
    gff = Gff(gff_file)

    gv = GenomeViz()
    track = gv.add_feature_track(gff.name, gff.range_size)
    track.add_gff_features(gff)

    fig_pngfile = tmp_path / "gff_feature.png"
    gv.savefig(fig_pngfile)
    assert fig_pngfile.exists()

    fig_htmlfile = fig_pngfile.with_suffix(".html")
    gv.savefig_html(fig_htmlfile)
    assert fig_htmlfile.exists()


def test_get_track_error():
    """Test get_track error"""
    gv = GenomeViz()
    with pytest.raises(ValueError) as e:
        gv.get_track("test")
    assert str(e.value).startswith("track_name='test' is not found")


def test_top_track_error():
    """Test top track property access error"""
    gv = GenomeViz()
    with pytest.raises(ValueError) as e:
        gv.top_track
    assert str(e.value).startswith("No track found.")


def test_max_track_size_error():
    """Test max track size property access error"""
    gv = GenomeViz()
    with pytest.raises(ValueError) as e:
        gv.max_track_size
    assert str(e.value).startswith("No track found.")


def test_track_name_duplication_error():
    """Test track name duplication case"""
    gv = GenomeViz()
    dup_name, size = "duplication track", 1000
    gv.add_feature_track(dup_name, size)
    with pytest.raises(ValueError):
        gv.add_feature_track(dup_name, size)


def test_add_link_for_no_adjacent_tracks():
    """Test add link to no adjacent feature tracks"""
    gv = GenomeViz()
    track_name1, track_name2, track_name3 = "track01", "track02", "track03"
    gv.add_feature_track(track_name1, 1000)
    gv.add_feature_track(track_name2, 1100)
    gv.add_feature_track(track_name3, 1200)

    gv.add_link((track_name1, 100, 200), (track_name2, 200, 300))
    with pytest.raises(ValueError):
        gv.add_link((track_name1, 100, 200), (track_name3, 200, 300))


def test_add_link_range_error():
    """Test add_link out range error"""
    gv = GenomeViz()
    track_name1, track_name2, track_name3 = "track01", "track02", "track03"
    gv.add_feature_track(track_name1, 1000)  # Range: 0 - 1000
    gv.add_feature_track(track_name2, 1100)  # Range: 0 - 1100
    gv.add_feature_track(track_name3, 1200, start_pos=100)  # Range: 100 - 1300

    # Case1. track_link1 is out range
    with pytest.raises(ValueError):
        gv.add_link((track_name1, 900, 1100), (track_name2, 200, 400))
    # Case2. track_link2 is out range
    with pytest.raises(ValueError):
        gv.add_link((track_name1, 800, 1000), (track_name2, 900, 1150))
    # Case3. track_link is out range (Set start_pos)
    with pytest.raises(ValueError):
        gv.add_link((track_name2, 800, 1000), (track_name3, 50, 250))
    with pytest.raises(ValueError):
        gv.add_link((track_name2, 800, 1000), (track_name3, 1150, 1350))


def test_top_track():
    """Test top track"""
    gv = GenomeViz()
    top_track_name, top_track_size = "Top Track", 1000
    gv.add_feature_track(top_track_name, top_track_size)
    gv.add_feature_track("track02", 1200)
    gv.add_feature_track("track03", 1100)

    assert gv.top_track.name == top_track_name
    assert gv.top_track.size == top_track_size


def test_bottom_track():
    """Test top track"""
    gv = GenomeViz()
    bottom_track_name, bottom_track_size = "Bottom Track", 1100
    gv.add_feature_track("track01", 1000)
    gv.add_feature_track("track02", 1200)
    gv.add_feature_track(bottom_track_name, bottom_track_size)

    assert gv.bottom_track.name == bottom_track_name
    assert gv.bottom_track.size == bottom_track_size


def test_ax_property():
    """Test ax property access"""
    gv = GenomeViz()
    gv.add_feature_track("track01", 1000)
    # Can't access ax property before calling 'plotfig' method
    with pytest.raises(ValueError):
        gv.top_track.ax
    # Propery access ax property after calling 'plotfig' method
    gv.plotfig()
    gv.top_track.ax


def test_offset_property():
    """Test offset property access"""
    gv = GenomeViz()
    gv.add_feature_track("track01", 1000)
    # Can't access offset property before calling 'plotfig' method
    with pytest.raises(ValueError):
        gv.top_track.offset
    # Propery access offset property after calling 'plotfig' method
    gv.plotfig()
    gv.top_track.offset
