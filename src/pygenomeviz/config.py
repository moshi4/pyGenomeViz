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


ann_adjust = _AnnotationAdjustConfig
