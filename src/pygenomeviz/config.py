from typing import ClassVar


class _AnnotationAdjustConfig:
    """Annotation Position Adjustment Config"""

    enabled: ClassVar[bool] = True
    """Enable annotation position adjustment (default: `True`)"""
    limit: ClassVar[int] = 200
    """Limit of track annotation number for position adjustment (default: `200`)"""
    dy: ClassVar[float] = 0.05
    """Delta y for iterative position adjustment (default: `0.05`)"""
    wpad: ClassVar[float] = 0.1
    """Text bbox width padded size for space between annotations (default: `0.1`)"""
    hpad: ClassVar[float] = 0.3
    """Text bbox height padded size for space between annotations (default: `0.3`)"""


class _PatchZorder:
    """Patch Zorder Config"""

    big_feature: ClassVar[float] = 2.0
    normal_feature: ClassVar[float] = 1.0
    intron: ClassVar[float] = 0.99
    link: ClassVar[float] = 1.0
    ann_line: ClassVar[float] = 0.99
    lollipop_line: ClassVar[float] = 0.99
    lollipop_point: ClassVar[float] = 2.01
    highlight: ClassVar[float] = 3.0

    @classmethod
    def feature(cls, bigstyle: bool) -> float:
        """Feature patch zorder"""
        return cls.big_feature if bigstyle else cls.normal_feature


ann_adjust = _AnnotationAdjustConfig
zorder = _PatchZorder
