import pytest
from pygenomeviz.feature import ExonFeature, Feature

###########################################################
# Test Feature Class
###########################################################


def test_feature_add():
    """Test feature __add__"""
    start, end = 10, 50
    feature = Feature(start, end, strand=1)
    add_value = 10
    new_feature = feature + add_value

    assert feature.start == start
    assert feature.end == end
    assert feature.length == end - start

    assert new_feature.start == start + add_value
    assert new_feature.end == end + add_value
    assert new_feature.length == end - start


def test_feature_check_error():
    """Test feature check error case"""
    # Case1. start > end
    with pytest.raises(ValueError):
        Feature(100, 0)
    # Case2. Invalid labelvpos
    with pytest.raises(ValueError):
        Feature(0, 100, labelvpos="invalid")  # type: ignore
    # Case3. Invalid labelhpos
    with pytest.raises(ValueError):
        Feature(0, 100, labelhpos="invalid")  # type: ignore
    # Case4. Invalid labelha
    with pytest.raises(ValueError):
        Feature(0, 100, labelha="invalid")  # type: ignore
    # Case5. Invalid plotstyle
    with pytest.raises(ValueError):
        Feature(0, 100, plotstyle="invalid")  # type: ignore
    # Case6. Invalid arrow_shaft_ratio range (not 0 <= v <= 1)
    with pytest.raises(ValueError):
        Feature(0, 100, arrow_shaft_ratio=2.0)
    # Case7. Invalid size_ratio range (not 0 <= v <= 1)
    with pytest.raises(ValueError):
        Feature(0, 100, size_ratio=2.0)


###########################################################
# Test Exon Feature Class
###########################################################


def test_exon_feature_add():
    """Test exon feature __add__"""
    exon_regions = [(10, 40), (50, 110), (130, 160)]
    exon_feature = ExonFeature(exon_regions)
    add_value = 10
    new_exon_feature = exon_feature + add_value

    assert exon_feature.start == 10
    assert exon_feature.end == 160
    assert exon_feature.length == 150

    assert new_exon_feature.start == exon_feature.start + add_value
    assert new_exon_feature.end == exon_feature.end + add_value
    assert new_exon_feature.length == 150
    assert new_exon_feature.exon_regions == [(20, 50), (60, 120), (140, 170)]


def test_exon_feature_intron_regions():
    """Test exon feature intron regions"""
    exon_regions = [(10, 40), (50, 110), (130, 160)]
    exon_feature = ExonFeature(exon_regions)
    expected_intron_regions = [(40, 50), (110, 130)]
    assert exon_feature.intron_regions == expected_intron_regions


def test_exon_feature_check_error():
    """Test exon feature check error case"""
    # Case1. exon_start > exon_end
    with pytest.raises(ValueError):
        exon_regions = [(40, 10), (110, 50), (160, 130)]
        ExonFeature(exon_regions)
    # Case2. exon regions is not ascending order of genomic regions
    with pytest.raises(ValueError):
        exon_regions = [(130, 160), (50, 110), (10, 40)]
        ExonFeature(exon_regions)
