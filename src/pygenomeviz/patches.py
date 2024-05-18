from __future__ import annotations

from matplotlib.patches import FancyArrow, PathPatch
from matplotlib.path import Path


class Arrow(FancyArrow):
    """Arrow Patch Class"""

    def __init__(
        self,
        start: int,
        end: int,
        strand: int,
        *,
        max_size: int,
        bigstyle: bool = False,
        show_head: bool = True,
        shaft_ratio: float = 0.5,
        ylim: tuple[float, float] = (-1, 1),
        **kwargs,
    ) -> None:
        """
        Parameters
        ----------
        start : int
            Start position
        end : int
            End position
        strand : int
            Strand
        max_size : int
            Axes x size (Required for calculation of arrow head length)
        bigstyle : bool, optional
            Big style or not
        shaft_ratio : float, optional
            Arrow shaft size ratio
        ylim : tuple[float, float], optional
            Axes y limit
        **kwargs : dict, optional
            Patch properties (e.g. `fc="red", ec="black", lw=1, hatch="//", ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        # Setup default patch kwargs
        default_zorder = 2 if bigstyle else 1
        kwargs.setdefault("zorder", default_zorder)
        kwargs.setdefault("lw", 0)
        kwargs.setdefault("clip_on", False)
        if "color" not in kwargs and "facecolor" not in kwargs:
            kwargs.setdefault("fc", "orange")

        # x, y
        x = end if strand == -1 else start
        if bigstyle:
            y = 0
        else:
            y = ylim[0] / 2 if strand == -1 else ylim[1] / 2
        # dx, dy
        length = end - start
        dx, dy = length * strand, 0
        # head width
        max_width = ylim[1] - ylim[0]
        head_width = max_width if bigstyle else max_width / 2
        # shaft_width
        shaft_width = head_width * shaft_ratio
        # head length
        head_length = max_size * 0.015
        if length < head_length:
            head_length = length

        if not show_head:
            head_length = 0
            head_width = shaft_width

        super().__init__(
            x,
            y,
            dx,
            dy,
            width=shaft_width,
            length_includes_head=True,
            head_width=head_width,
            head_length=head_length,
            **kwargs,
        )


class BigArrow(Arrow):
    """BigArrow Patch Class"""

    def __init__(
        self,
        start: int,
        end: int,
        strand: int,
        *,
        max_size: int,
        show_head: bool = True,
        shaft_ratio: float = 0.5,
        ylim: tuple[float, float] = (-1, 1),
        **kwargs,
    ) -> None:
        super().__init__(
            start,
            end,
            strand,
            max_size=max_size,
            bigstyle=True,
            show_head=show_head,
            shaft_ratio=shaft_ratio,
            ylim=ylim,
            **kwargs,
        )


class Box(PathPatch):
    """Box Patch Class"""

    def __init__(
        self,
        start: int,
        end: int,
        strand: int,
        *,
        bigstyle: bool = False,
        ylim: tuple[float, float] = (-1, 1),
        **kwargs,
    ) -> None:
        """
        Parameters
        ----------
        start : int
            Start position
        end : int
            End position
        strand : int
            Strand
        bigstyle : bool, optional
            Big style or not
        ylim : tuple[float, float], optional
            Axes y limit
        **kwargs : dict, optional
            Patch properties (e.g. `fc="red", ec="black", lw=1, hatch="//", ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        # Setup default patch kwargs
        default_zorder = 2 if bigstyle else 1
        kwargs.setdefault("zorder", default_zorder)
        kwargs.setdefault("lw", 0)
        kwargs.setdefault("clip_on", False)
        if "color" not in kwargs and "facecolor" not in kwargs:
            kwargs.setdefault("fc", "orange")

        if bigstyle:
            lower_y, upper_y = ylim[0], ylim[1]
        else:
            lower_y, upper_y = (ylim[0], 0) if strand == -1 else (0, ylim[1])

        p1 = (start, lower_y)
        p2 = (end, lower_y)
        p3 = (end, upper_y)
        p4 = (start, upper_y)
        box = Path([p1, p2, p3, p4, p1], closed=True)  # type: ignore

        super().__init__(box, **kwargs)


class BigBox(Box):
    """BigBox Patch Class"""

    def __init__(
        self,
        start: int,
        end: int,
        strand: int,
        *,
        ylim: tuple[float, float] = (-1, 1),
        **kwargs,
    ) -> None:
        super().__init__(
            start,
            end,
            strand,
            bigstyle=True,
            ylim=ylim,
            **kwargs,
        )


class RoundBox(PathPatch):
    """RoundBox Patch Class"""

    def __init__(
        self,
        start: int,
        end: int,
        strand: int,
        *,
        max_size: int,
        bigstyle: bool = False,
        ylim: tuple[float, float] = (-1, 1),
        **kwargs,
    ) -> None:
        """
        Parameters
        ----------
        start : int
            Start position
        end : int
            End position
        strand : int
            Strand
        max_size : int
            Axes x size (Required for calculation of round box radius)
        bigstyle : bool, optional
            Big style or not
        ylim : tuple[float, float], optional
            Axes y limit
        **kwargs : dict, optional
            Patch properties (e.g. `fc="red", ec="black", lw=1, hatch="//", ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        # Set default patch properties
        default_zorder = 2 if bigstyle else 1
        kwargs.setdefault("zorder", default_zorder)
        kwargs.setdefault("lw", 0)
        kwargs.setdefault("clip_on", False)
        if "color" not in kwargs and "facecolor" not in kwargs:
            kwargs.setdefault("fc", "orange")

        length = end - start
        r_size = max_size * 0.005
        if length <= 4 * r_size:
            r_size = length * 0.25

        xmin, xmax = start, end
        if bigstyle:
            ymin, ymax, ycenter = ylim[0], ylim[1], 0
        else:
            if strand == -1:
                ymin, ymax, ycenter = ylim[0], 0, ylim[0] / 2
            else:
                ymin, ymax, ycenter = 0, ylim[1], ylim[1] / 2

        path_data = [
            (Path.MOVETO, (xmin + r_size, ymax)),
            (Path.LINETO, (xmax - r_size, ymax)),
            (Path.CURVE3, (xmax + r_size, ycenter)),
            (Path.CURVE3, (xmax - r_size, ymin)),
            (Path.LINETO, (xmin + r_size, ymin)),
            (Path.CURVE3, (xmin - r_size, ycenter)),
            (Path.CURVE3, (xmin + r_size, ymax)),
        ]
        codes, verts = zip(*path_data)
        super().__init__(Path(verts, codes), **kwargs)


class BigRoundBox(RoundBox):
    """BigRoundBox Patch Class"""

    def __init__(
        self,
        start: int,
        end: int,
        strand: int,
        *,
        max_size: int,
        ylim: tuple[float, float] = (-1, 1),
        **kwargs,
    ) -> None:
        super().__init__(
            start,
            end,
            strand,
            bigstyle=True,
            max_size=max_size,
            ylim=ylim,
            **kwargs,
        )


class Intron(PathPatch):
    """Intron Patch Class"""

    def __init__(
        self,
        start: int,
        end: int,
        strand: int,
        *,
        bigstyle: bool = False,
        ylim: tuple[float, float] = (-1, 1),
        **kwargs,
    ) -> None:
        """
        Parameters
        ----------
        start : int
            Start position
        end : int
            End position
        strand : int
            Strand
        bigstyle : bool, optional
            Big style or not
        ylim : tuple[float, float], optional
            Axes y limit
        **kwargs : dict, optional
            Patch properties (e.g. `fc="red", ec="black", lw=1, hatch="//", ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        # Setup default patch kwargs
        kwargs.setdefault("zorder", 0.99)
        kwargs.setdefault("lw", 1)
        kwargs.setdefault("clip_on", False)
        kwargs.setdefault("color", "black")
        kwargs.update(fill=False)

        xmin, xmax = min(start, end), max(start, end)
        xcenter = (xmin + xmax) / 2
        if bigstyle:
            ymin, ymax, ycenter = ylim[0], ylim[1], 0
        else:
            if strand == -1:
                ymin, ymax, ycenter = ylim[0], 0, ylim[0] / 2
            else:
                ymin, ymax, ycenter = 0, ylim[1], ylim[1] / 2
        ytop = ymin if strand == -1 else ymax

        path_data = [
            (Path.MOVETO, (xmin, ycenter)),
            (Path.LINETO, (xcenter, ytop)),
            (Path.LINETO, (xmax, ycenter)),
        ]
        codes, verts = zip(*path_data)
        super().__init__(Path(verts, codes), **kwargs)


class Link(PathPatch):
    """Link Patch Class"""

    def __init__(
        self,
        start1: int,
        end1: int,
        start2: int,
        end2: int,
        *,
        ylim: tuple[float, float] = (-1, 1),
        curve: bool = False,
        **kwargs,
    ):
        """
        Parameters
        ----------
        start1 : int
            Link upper start postion
        end1 : int
            Link upper end position
        start2 : int
            Link lower start position
        end2 : int
            Link lower end position
        ylim : tuple[float, float], optional
            Axes y limit
        curve : bool, optional
            Curve or not
        **kwargs : dict, optional
            Patch properties (e.g. `fc="red", ec="black", lw=1, hatch="//", ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        # Set default patch properties
        kwargs.setdefault("lw", 0.1)  # If lw=0, twisted part is almost invisible
        kwargs.setdefault("zorder", 1)
        kwargs.setdefault("clip_on", False)
        if "color" not in kwargs and "facecolor" not in kwargs:
            kwargs.setdefault("fc", "grey")

        ymin, ymax = min(ylim), max(ylim)

        if curve:
            ctl_y_point1, ctl_y_point2 = ymax / 3, ymin / 3
            path_data = [
                (Path.MOVETO, (start2, ymin)),
                (Path.LINETO, (end2, ymin)),
                (Path.CURVE4, (end2, ctl_y_point2)),
                (Path.CURVE4, (end1, ctl_y_point1)),
                (Path.LINETO, (end1, ymax)),
                (Path.LINETO, (start1, ymax)),
                (Path.CURVE4, (start1, ctl_y_point1)),
                (Path.CURVE4, (start2, ctl_y_point2)),
                (Path.LINETO, (start2, ymin)),
            ]
        else:
            path_data = [
                (Path.MOVETO, (start2, ymin)),
                (Path.LINETO, (end2, ymin)),
                (Path.LINETO, (end1, ymax)),
                (Path.LINETO, (start1, ymax)),
                (Path.LINETO, (start2, ymin)),
            ]
        codes, verts = zip(*path_data)
        super().__init__(Path(verts, codes), **kwargs)


PLOTSTYLE2PATCH = dict(
    arrow=Arrow,
    bigarrow=BigArrow,
    box=Box,
    bigbox=BigBox,
    rbox=RoundBox,
    bigrbox=BigRoundBox,
)
