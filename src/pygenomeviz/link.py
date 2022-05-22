from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple


@dataclass
class Link:
    track_name1: str
    track_start1: int
    track_end1: int
    track_name2: str
    track_start2: int
    track_end2: int
    normal_color: str = "grey"
    inverted_color: str = "red"

    @property
    def color(self) -> str:
        """Get color"""
        if self.is_inverted():
            return self.inverted_color
        else:
            return self.normal_color

    def get_link_pos(self, track_name: str) -> Tuple[int, int]:
        """Get link start-end position tuple"""
        if self.track_name1 == track_name:
            return (self.track_start1, self.track_end1)
        elif self.track_name2 == track_name:
            return (self.track_start2, self.track_end2)
        else:
            raise ValueError()

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
        )
