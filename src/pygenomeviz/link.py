from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from typing import Any, Dict, Optional, Tuple

from matplotlib import colors
from matplotlib.figure import Axes
from matplotlib.patches import PathPatch
from matplotlib.path import Path


@dataclass
class Link:
    """Link DataClass"""

    track_name1: str
    track_start1: int
    track_end1: int
    track_name2: str
    track_start2: int
    track_end2: int
    normal_color: str = "grey"
    inverted_color: str = "red"
    alpha: float = 0.8
    v: Optional[float] = None
    vmin: float = 0
    vmax: float = 100
    curve: bool = False
    size_ratio: float = 0.9
    patch_kws: Optional[Dict[str, Any]] = None

    def __post_init__(self):
        # Check color string
        if not colors.is_color_like(self.normal_color):
            err_msg = f"'normal_color={self.normal_color}' is not color like."
            raise ValueError(err_msg)
        if not colors.is_color_like(self.inverted_color):
            err_msg = f"'inverted_color={self.inverted_color}' is not color like."
            raise ValueError(err_msg)
        # Check size ratio
        if not 0 <= self.size_ratio <= 1:
            err_msg = "'size_ratio' must be '0 <= value <= 1' "
            err_msg += f"(value={self.size_ratio})"
            raise ValueError(err_msg)
        # Check range of interpolation value
        if self.v is not None:
            if not 0 <= self.vmin <= self.v <= self.vmax <= 100:
                err_msg = "'Interpolation value must be "
                err_msg += "'0 <= vmin <= value <= vmax <= 100' "
                err_msg += f" (value={self.v}, vmin={self.vmin}, vmax={self.vmax})"
                raise ValueError(err_msg)

    def plot_link(self, ax: Axes, ylim: Tuple[float, float] = (-1.0, 1.0)) -> None:
        """Plot link

        Parameters
        ----------
        ax : Axes
            Matplotlib axes object to be plotted
        ylim: Tuple[flaot, float], optional
            Min-Max y coordinatess
        """
        ymin, ymax = ylim[0] * self.size_ratio, ylim[1] * self.size_ratio
        start1, end1 = self.track_start1 - 1, self.track_end1
        start2, end2 = self.track_start2 - 1, self.track_end2
        if self.curve:
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
        ax.add_patch(PathPatch(Path(verts, codes), **self._patch_kwargs()))

    def _patch_kwargs(self) -> Dict[str, Any]:
        """Patch keyword arguments dict

        Returns
        -------
        patch_kwargs : Dict[str, Any]
            Patch keyword arguments dict
        """
        patch_kws = {} if self.patch_kws is None else self.patch_kws
        return {
            "fc": self.color,
            "lw": 0,
            **patch_kws,
        }

    @property
    def color(self) -> str:
        """Get conditional hexcolor code"""
        color = self.inverted_color if self.is_inverted else self.normal_color
        if self.v is None:
            rgba = colors.to_rgba(color, alpha=self.alpha)
            return colors.to_hex(rgba, keep_alpha=True)
        else:

            def to_nearly_white(color: str, nearly_value: float = 0.1) -> str:
                """Convert target color to nearly white"""
                cmap = colors.LinearSegmentedColormap.from_list("m", ("white", color))
                return colors.to_hex(cmap(nearly_value))

            nearly_white = to_nearly_white(color)
            cmap = colors.LinearSegmentedColormap.from_list("m", (nearly_white, color))
            norm = colors.Normalize(vmin=self.vmin, vmax=self.vmax)
            norm_value = norm(self.v)
            return colors.to_hex(cmap(norm_value, alpha=self.alpha), keep_alpha=True)

    @property
    def is_inverted(self) -> bool:
        """Check inverted link or not

        Returns
        -------
        is_inverted : bool
            Check result
        """
        track_link_length1 = self.track_end1 - self.track_start1
        track_link_length2 = self.track_end2 - self.track_start2
        return track_link_length1 * track_link_length2 < 0

    def add_offset(self, track_name2offset: Dict[str, int]) -> Link:
        """Add offset to each link position

        Parameters
        ----------
        track_name2offset : Dict[str, int]
            Track name & offset dict

        Returns
        -------
        link : Link
            Offset added Link object
        """
        link = deepcopy(self)
        link.track_start1 += track_name2offset[self.track_name1]
        link.track_end1 += track_name2offset[self.track_name1]
        link.track_start2 += track_name2offset[self.track_name2]
        link.track_end2 += track_name2offset[self.track_name2]
        return link
