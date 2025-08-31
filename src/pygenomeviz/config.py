from typing import ClassVar


class _AnnotationAdjustConfig:
    """Annotation Position Adjustment Config"""

    enabled: ClassVar[bool] = True
    """Enable annotation position adjustment (default: `True`)"""
    limit: ClassVar[int] = 500
    """Limit of track annotation number for position adjustment (default: `500`)"""
    dy: ClassVar[float] = 0.1
    """Delta y for iterative position adjustment (default: `0.1`)"""
    expand: ClassVar[tuple[float, float]] = (1.2, 1.5)
    """Expand width & height factor of text bbox (default: `(1.2, 1.5)`)"""


ann_adjust = _AnnotationAdjustConfig
