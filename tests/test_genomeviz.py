from __future__ import annotations

from pathlib import Path

import pytest
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation

from pygenomeviz import GenomeViz
from pygenomeviz.parser import Genbank, Gff


def test_plot_features_with_all_plotstyle(tmp_path: Path):
    """Test plot features with all plotstyle"""
    gv = GenomeViz()
    gv.set_scale_bar()
    gv.set_scale_xticks()
    track = gv.add_feature_track("test", 1000)

    track.add_feature(start=50, end=150, strand=1)
    track.add_feature(start=200, end=300, strand=-1, plotstyle="arrow")
    track.add_feature(start=330, end=400, strand=1, plotstyle="bigbox")
    track.add_feature(start=420, end=500, strand=1, plotstyle="box")
    track.add_feature(start=550, end=600, strand=1, plotstyle="bigrbox")
    track.add_feature(start=650, end=750, strand=-1, plotstyle="rbox")
    track.add_feature(start=780, end=880, strand=1, plotstyle="bigarrow")

    gv.savefig(tmp_path / "result.png")
    gv.savefig_html(tmp_path / "result.html")


def test_exon_features_plot(tmp_path: Path):
    """Test exon features plot"""
    gv = GenomeViz()
    track = gv.add_feature_track("track1", 1000)
    locs = [(0, 100), (200, 300), (350, 400)]
    track.add_exon_feature(locs)

    feature = SeqFeature(
        location=CompoundLocation(
            [
                SimpleLocation(500, 550, 1),
                SimpleLocation(600, 700, 1),
                SimpleLocation(800, 900, 1),
            ]
        )
    )
    track.add_exon_features(feature, plotstyle="arrow")

    gv.savefig(tmp_path / "result.png")
    gv.savefig_html(tmp_path / "result.html")


def test_manual_features_links_plot(tmp_path: Path):
    """Test manual features & links plot"""
    genome_list = [
        dict(
            name="genome 01",
            size=1000,
            features=((150, 300, 1), (500, 700, -1), (750, 950, 1)),
        ),
        dict(
            name="genome 02",
            size=1300,
            features=((50, 200, 1), (350, 450, 1), (700, 900, -1), (950, 1150, -1)),
        ),
        dict(
            name="genome 03",
            size=1200,
            features=((150, 300, 1), (350, 450, -1), (500, 700, -1), (700, 900, -1)),
        ),
    ]

    gv = GenomeViz(track_align_type="center")
    for genome in genome_list:
        name, size, features = genome["name"], genome["size"], genome["features"]
        track = gv.add_feature_track(name, size, align_label=True)
        for idx, feature in enumerate(features, 1):
            start, end, strand = feature
            track.add_feature(start, end, strand, label=f"gene{idx:02d}", lw=1)

    # Add links between "genome 01" and "genome 02"
    gv.add_link(("genome 01", 150, 300), ("genome 02", 50, 200))
    gv.add_link(("genome 01", 700, 500), ("genome 02", 900, 700))
    gv.add_link(("genome 01", 750, 950), ("genome 02", 1150, 950))
    # Add links between "genome 02" and "genome 03"
    kwargs = dict(color="skyblue", inverted_color="lime", curve=True)
    gv.add_link(("genome 02", 50, 200), ("genome 03", 150, 300), **kwargs)  # type: ignore
    gv.add_link(("genome 02", 350, 450), ("genome 03", 450, 350), **kwargs)  # type: ignore
    gv.add_link(("genome 02", 900, 700), ("genome 03", 700, 500), **kwargs)  # type: ignore
    gv.add_link(("genome 03", 900, 700), ("genome 02", 1150, 950), **kwargs)  # type: ignore

    gv.savefig(tmp_path / "result.png")
    gv.savefig_html(tmp_path / "result.html")


def test_genbank_plot(gbk_file: Path, tmp_path: Path):
    """Test genbank features plot"""
    gbk = Genbank(gbk_file)

    gv = GenomeViz()
    gv.set_scale_xticks()
    track = gv.add_feature_track(gbk.name, gbk.genome_length)
    features = gbk.extract_features()
    track.add_features(features)

    gv.savefig(tmp_path / "result.png")
    gv.savefig_html(tmp_path / "result.html")


def test_gff_plot(gff_file: Path, tmp_path: Path):
    """Test gff features plot"""
    gff = Gff(gff_file)

    gv = GenomeViz()
    gv.set_scale_xticks()
    track = gv.add_feature_track(gff.name, gff.genome_length)
    features = gff.extract_features()
    track.add_features(features)

    gv.savefig(tmp_path / "result.png")
    gv.savefig_html(tmp_path / "result.html")


def test_savefig_html_failed(gbk_file: Path, tmp_path: Path):
    """Test `gv.savefig_html()` failed when fast_render=True"""
    gbk = Genbank(gbk_file)

    gv = GenomeViz()
    track = gv.add_feature_track(gbk.name, gbk.genome_length)
    track.add_features(gbk.extract_features())
    fig = gv.plotfig(fast_render=True)
    with pytest.raises(ValueError):
        gv.savefig_html(tmp_path / "result.html", fig)
