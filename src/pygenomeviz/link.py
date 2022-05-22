from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

from matplotlib.colors import LinearSegmentedColormap, Normalize, to_hex


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
    identity: Optional[float] = None
    interpolation: bool = True

    @property
    def color(self) -> str:
        """Get color"""
        color = self.inverted_color if self.is_inverted() else self.normal_color
        if self.interpolation is False or self.identity is None:
            return color
        else:
            cmap = LinearSegmentedColormap.from_list("cmap", ("white", color))
            norm = Normalize(vmin=0, vmax=100)
            norm_value = norm(self.identity)
            return to_hex(cmap(norm_value))

    def is_inverted(self):
        """Check inverted link or not"""
        track_link_length1 = self.track_end1 - self.track_start1
        track_link_length2 = self.track_end2 - self.track_start2
        return track_link_length1 * track_link_length2 < 0

    def add_offset(self, track_name2offset: Dict[str, int]) -> Link:
        """Add offset to each link position"""
        return Link(
            self.track_name1,
            self.track_start1 + track_name2offset[self.track_name1],
            self.track_end1 + track_name2offset[self.track_name1],
            self.track_name2,
            self.track_start2 + track_name2offset[self.track_name2],
            self.track_end2 + track_name2offset[self.track_name2],
            self.normal_color,
            self.inverted_color,
            self.identity,
            self.interpolation,
        )
