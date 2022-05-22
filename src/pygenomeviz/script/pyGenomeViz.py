from pygenomeviz.genomeviz import GenomeViz


def main():
    """Test main script"""
    gv = GenomeViz(fig_width=15, fig_track_height=1.0, align_type="center")
    # Track01
    track1 = gv.add_track(name="track1", size=100)
    track1.add_feature(20, 30, 1)
    track1.add_feature(50, 55, -1)
    # Track02
    track2 = gv.add_track(name="track2", size=120)
    track2.add_feature(20, 40, 1, facecolor="blue", edgecolor="blue")
    track2.add_feature(50, 80, -1, plotstyle="box", edgecolor="black")
    # Track03
    track3 = gv.add_track(name="track3", size=200)
    track3.add_feature(10, 15, -1)
    track3.add_feature(20, 70, 1, plotstyle="arrow")
    track3.add_feature(100, 102, -1, plotstyle="arrow")
    track3.add_feature(135, 136, 1, plotstyle="box")
    # Track04
    track4 = gv.add_track(name="track4", size=170)
    track4.add_feature(10, 15, -1)
    track4.add_feature(20, 70, 1, plotstyle="arrow")
    track4.add_feature(100, 102, -1, plotstyle="arrow")
    track4.add_feature(135, 136, 1, plotstyle="box", edgecolor="black")

    # Track01 - Track02 link
    gv.add_link(("track1", 10, 20), ("track2", 0, 10))
    gv.add_link(("track2", 50, 70), ("track3", 80, 110))
    gv.add_link(("track3", 50, 20), ("track4", 0, 30))

    gv.plotfig("example/test.jpg")


if __name__ == "__main__":
    main()
