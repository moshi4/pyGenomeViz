from __future__ import annotations

import pytest

from pygenomeviz.exception import FeatureRangeError
from pygenomeviz.track import FeatureTrack


def test_feature_track_properties():
    """Test FeatureTrack properties"""
    track_name = "feature_track"
    seg_name2range = dict(track1=(0, 100), track2=(200, 1000), track3=(2000, 2500))
    track = FeatureTrack(track_name, seg_name2range)

    assert track.name == track_name
    for seg, (seg_name, range) in zip(track.segments, seg_name2range.items()):
        assert seg.name == seg_name
        assert (seg.start, seg.end) == range


def test_add_outside_range_feature():
    """Test add outside range feature"""
    track = FeatureTrack("test", dict(a=(0, 1000), b=(1000, 2000)))
    # Success
    track.add_feature(200, 800, 1)
    # Raise outside range error in first segment
    with pytest.raises(FeatureRangeError):
        track.add_feature(800, 1500, -1)
    # Raise outside range error in second segment
    with pytest.raises(FeatureRangeError):
        track.add_exon_feature([(500, 1500)], -1, target_seg="b")
