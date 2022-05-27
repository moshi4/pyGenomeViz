from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

from matplotlib.colors import LinearSegmentedColormap, Normalize, is_color_like, to_hex


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
    interpolation_value: Optional[float] = None
    vmin: float = 0
    vmax: float = 100

    def __post_init__(self):
        if not is_color_like(self.normal_color):
            err_msg = f"'normal_color={self.normal_color}' is not color like."
            raise ValueError(err_msg)
        if not is_color_like(self.inverted_color):
            err_msg = f"'inverted_color={self.inverted_color}' is not color like."
            raise ValueError(err_msg)
        if self.interpolation_value is not None:
            if not self.vmin <= self.interpolation_value <= self.vmax:
                err_msg = "'Interpolation value must be "
                err_msg += f"'{self.vmin} <= value <= {self.vmax}' "
                err_msg += f"(value={self.interpolation_value})"
                raise ValueError(err_msg)

    def polygon_xy(
        self, ymin: float = -1.0, ymax: float = 1.0
    ) -> Tuple[
        Tuple[float, float],
        Tuple[float, float],
        Tuple[float, float],
        Tuple[float, float],
    ]:
        """Plygon coordinate for link plot

        Parameters
        ----------
        ymin : float, optional
            Polygon min y coordinate
        ymax : float, optional
            Polygon max y coordinate

        Returns
        -------
        polygon_xy : tuple
            Polygon xy coordinates
        """
        return (
            (self.track_start2, ymin),
            (self.track_end2, ymin),
            (self.track_end1, ymax),
            (self.track_start1, ymax),
        )

    @property
    def color(self) -> str:
        """Get conditional hexcolor code

        If `self.interpolation_value` is set, interpolated color is calculated.
        """
        color = self.inverted_color if self.is_inverted() else self.normal_color
        if self.interpolation_value is None:
            return color
        else:
            cmap = LinearSegmentedColormap.from_list("cmap", ("white", color))
            norm = Normalize(vmin=self.vmin, vmax=self.vmax)
            norm_value = norm(self.interpolation_value)
            return to_hex(cmap(norm_value))

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
        return Link(
            self.track_name1,
            self.track_start1 + track_name2offset[self.track_name1],
            self.track_end1 + track_name2offset[self.track_name1],
            self.track_name2,
            self.track_start2 + track_name2offset[self.track_name2],
            self.track_end2 + track_name2offset[self.track_name2],
            self.normal_color,
            self.inverted_color,
            self.interpolation_value,
            self.vmin,
            self.vmax,
        )
