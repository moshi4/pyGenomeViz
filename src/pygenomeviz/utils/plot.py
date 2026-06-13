from __future__ import annotations

from collections import defaultdict
from typing import TYPE_CHECKING

from matplotlib.collections import PatchCollection

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.patches import Patch


def plot_patches(
    patches: list[Patch],
    ax: Axes,
    fast_render: bool = True,
) -> None:
    """Plot patches

    Parameters
    ----------
    patches : list[Patch]
        Patches
    ax : Axes
        Axes
    fast_render: bool, optional
        Enable fast rendering using PatchCollection plot style.
    """
    if fast_render:
        zorder2patches: dict[float, list[Patch]] = defaultdict(list)
        for patch in patches:
            if patch.get_hatch() is not None:
                ax.add_patch(patch)
            else:
                zorder2patches[float(patch.get_zorder())].append(patch)
        for zorder, patch_group in zorder2patches.items():
            patch_col = PatchCollection(
                patch_group, match_original=True, zorder=zorder, clip_on=False
            )
            ax.add_collection(patch_col)
    else:
        for patch in patches:
            ax.add_patch(patch)
