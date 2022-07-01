import pytest
from pygenomeviz.track import FeatureTrack


def test_add_subtrack():
    """Test add subtrack"""
    track_name, track_size = "feature track", 1000
    track = FeatureTrack(track_name, track_size)
    track.add_subtrack()
    track.add_subtrack(ratio=0.7)

    assert len(track.subtracks) == 2
    assert track.subtracks[0].name == f"{track_name}_subtrack1"
    assert track.subtracks[0].size == track_size
    assert track.subtracks[0].ratio == 1.0
    assert track.subtracks[1].name == f"{track_name}_subtrack2"
    assert track.subtracks[1].size == track_size
    assert track.subtracks[1].ratio == 0.7


def test_add_feature_range_error():
    """Test add feature range check"""
    track_name, track_size = "test", 1000
    under_track_size, over_track_size = -100, track_size + 100
    track = FeatureTrack(track_name, track_size)
    # Case1. start < 0
    with pytest.raises(ValueError):
        track.add_feature(under_track_size, 100)
    # Case2: end > track_size
    with pytest.raises(ValueError):
        track.add_feature(100, over_track_size)
    # Case3: end > start
    with pytest.raises(ValueError):
        track.add_feature(500, 400)
    # Case4: start < 0 and end > track_size
    with pytest.raises(ValueError):
        track.add_feature(under_track_size, over_track_size)


def test_add_exon_feature_range_error():
    """Test add exon feature range check"""
    track_name, track_size = "test", 1000
    under_track_size, over_track_size = -100, track_size + 100
    track = FeatureTrack(track_name, track_size)
    # Case1. start < 0
    with pytest.raises(ValueError):
        exon_regions = [(under_track_size, 200), (300, 600), (650, 900)]
        track.add_exon_feature(exon_regions)
    # Case2: end > track_size
    with pytest.raises(ValueError):
        exon_regions = [(100, 200), (300, 600), (650, over_track_size)]
        track.add_exon_feature(exon_regions)
    # Case3: end > start
    with pytest.raises(ValueError):
        exon_regions = [(100, 200), (600, 300), (650, 900)]
        track.add_exon_feature(exon_regions)
    # Case4: start < 0 and end > track_size
    with pytest.raises(ValueError):
        exon_regions = [(under_track_size, 200), (600, 300), (650, over_track_size)]
        track.add_exon_feature(exon_regions)
