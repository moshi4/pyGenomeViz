from __future__ import annotations

from Bio.SeqFeature import SeqFeature, SimpleLocation

from pygenomeviz.utils import is_pseudo_feature, to_stack_features


def test_is_pseudo_feature():
    """Test `is_pseudo_feature()`"""
    no_pseudo_feature = SeqFeature(SimpleLocation(0, 100))
    assert is_pseudo_feature(no_pseudo_feature) is False
    pseudo_feature = SeqFeature(SimpleLocation(0, 100), qualifiers=dict(pseudo=[]))
    assert is_pseudo_feature(pseudo_feature) is True


def test_to_stack_features():
    """Test `to_stack_features()`"""
    locs = [(10, 20), (15, 25), (18, 30), (30, 50), (70, 100), (80, 85), (90, 100)]
    features = [SeqFeature(SimpleLocation(loc[0], loc[1])) for loc in locs]

    def to_loc(feature: SeqFeature) -> tuple[int, int]:
        """Convert feature to (start, end) location"""
        return (int(feature.location.start), int(feature.location.end))  # type: ignore

    stack_locs = []
    for sublist_features in to_stack_features(features):
        stack_locs.append(list(map(to_loc, sublist_features)))

    assert stack_locs == [
        [(10, 20), (30, 50), (70, 100)],
        [(15, 25), (80, 85), (90, 100)],
        [(18, 30)],
    ]
