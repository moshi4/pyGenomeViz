from __future__ import annotations

from typing import overload

import matplotlib.pyplot as plt
import numpy as np
from Bio.SeqFeature import SeqFeature
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_hex

from pygenomeviz.typing import Unit


class ColorCycler:
    """Color Cycler Class"""

    counter = 0
    cmap = plt.get_cmap("tab10")  # type: ignore

    def __new__(cls, n: int | None = None) -> str:
        """Get hexcolor cyclically from cmap by counter or user specified number

        `ColorCycler()` works same as `ColorCycler.get_color()` (syntactic sugar)

        Parameters
        ----------
        n : int | None, optional
            Number for color cycle. If None, counter class variable is used.

        Returns
        -------
        hexcolor : str
            Cyclic hexcolor string
        """
        return cls.get_color(n)

    @classmethod
    def reset_cycle(cls) -> None:
        """Reset cycle counter"""
        cls.counter = 0

    @classmethod
    def set_cmap(cls, name: str) -> None:
        """Set colormap (Default: `tab10`)"""
        cls.cmap = plt.get_cmap(name)  # type: ignore
        cls.counter = 0

    @classmethod
    def get_color(cls, n: int | None = None) -> str:
        """Get hexcolor cyclically from cmap by counter or user specified number

        Parameters
        ----------
        n : int | None, optional
            Number for color cycle. If None, counter class variable is used.

        Returns
        -------
        hexcolor : str
            Cyclic hexcolor string
        """
        if n is None:
            n = cls.counter
            cls.counter += 1
        return to_hex(cls.cmap(n % cls.cmap.N), keep_alpha=True)

    @classmethod
    def get_color_list(cls, n: int | None = None) -> list[str]:
        """Get hexcolor list of colormap

        Parameters
        ----------
        n : int | None, optional
            If n is None, all(=cmap.N) hexcolors are extracted from colormap.
            If n is specified, hexcolors are extracted from n equally divided colormap.

        Returns
        -------
        hexcolor_list : list[str]
            Hexcolor list
        """
        if n is None:
            cmap_idx_list = list(range(0, cls.cmap.N))
        elif n > 0:
            cmap_idx_list = [int(i) for i in np.linspace(0, cls.cmap.N, n)]
        else:
            raise ValueError(f"{n=} is invalid number (Must be 'n > 0').")

        return [to_hex(cls.cmap(i), keep_alpha=True) for i in cmap_idx_list]


def is_pseudo_feature(feature: SeqFeature) -> bool:
    """Check target feature is pseudo or not from qualifiers tag

    Parameters
    ----------
    feature : SeqFeature
        Target feature

    Returns
    -------
    check_result : bool
        pseudo check result
    """
    quals = feature.qualifiers
    return True if "pseudo" in quals or "pseudogene" in quals else False


def extract_features_within_range(
    features: list[SeqFeature],
    *,
    target_range: tuple[int, int],
) -> list[SeqFeature]:
    """Extract features by target range

    Parameters
    ----------
    features : list[SeqFeature]
        Features to be extracted
    target_range : tuple[int, int]
        Target range

    Returns
    -------
    range_features : list[SeqFeature]
        Features within target range
    """
    range_features: list[SeqFeature] = []
    for feature in features:
        start = int(feature.location.start)  # type: ignore
        end = int(feature.location.end)  # type: ignore
        if min(target_range) <= start <= end <= max(target_range):
            range_features.append(feature)
    return range_features


def interpolate_color(
    base_color: str,
    v: float,
    vmin: float = 0,
    vmax: float = 100,
) -> str:
    """Interpolate the base color between vmin and vmax

    `vmin[nearly white] <= v <= vmax[base_color]`

    Parameters
    ----------
    base_color : str
        Base color for interpolation
    v : float
        Interpolation value
    vmin : float, optional
        Min value
    vmax : float, optional
        Max value

    Returns
    -------
    interpolate_color : str
        Interpolated hexcolor
    """
    if not vmin <= v <= vmax:
        raise ValueError(f"{v=} is out of range ({vmin=}, {vmax=})")

    def to_nearly_white(color: str, nearly_value: float = 0.1) -> str:
        """Convert target color to nearly white"""
        cmap = LinearSegmentedColormap.from_list("m", ["white", color])
        return to_hex(cmap(nearly_value))  # type: ignore

    nearly_white = to_nearly_white(base_color)
    cmap = LinearSegmentedColormap.from_list("m", [nearly_white, base_color])
    norm = Normalize(vmin=vmin, vmax=vmax)
    return to_hex(cmap(norm(v)))  # type: ignore


@overload
def size_label_formatter(size: float, unit: Unit | None = None) -> str: ...
@overload
def size_label_formatter(size: list[float], unit: Unit | None = None) -> list[str]: ...
def size_label_formatter(
    size: float | list[float], unit: Unit | None = None
) -> str | list[str]:
    """Format scale size to human readable style (e.g. 1000 -> `1.0 Kb`))

    Parameters
    ----------
    size : float | list[float]
        Scale size (or size list)
    unit : Unit | None, optional
        Format target unit (`Gb`|`Mb`|`Kb`|`bp`)

    Returns
    -------
    format_size : str | list[str]
        Formatted size (or size list)
    """
    # Convert to list value
    if not isinstance(size, list):
        size = [size]

    # Get target unit & unit size
    # unit2unit_size: dict[Unit, int] = dict(Gb=10**9, Mb=10**6, Kb=10**3, bp=1)
    unit2unit_size: dict[Unit, int] = {"Gb": 10**9, "Mb": 10**6, "Kb": 10**3, "bp": 1}
    if unit is None:
        max_size = max(size)
        target_unit, target_unit_size = "", 0
        for unit, unit_size in unit2unit_size.items():
            if max_size >= unit_size:
                target_unit, target_unit_size = unit, unit_size
                break
    else:
        if unit not in unit2unit_size:
            raise ValueError(f"{unit=} is invalid (Gb, Mb, Kb, bp)")
        else:
            target_unit, target_unit_size = unit, unit2unit_size[unit]

    # Format size
    format_size: list[str] = []
    for value in size:
        if value == 0:
            format_size.append(f"0 {target_unit}")
        else:
            plot_size = value / target_unit_size
            str_format = ",.0f" if plot_size >= 10 else ".1f"
            format_size.append(f"{value / target_unit_size:{str_format}} {target_unit}")

    return format_size[0] if len(format_size) == 1 else format_size
