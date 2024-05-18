from __future__ import annotations

from collections import defaultdict

from matplotlib.axes import Axes
from matplotlib.collections import PatchCollection
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
        for zorder, patches in zorder2patches.items():
            patch_col = PatchCollection(
                patches, match_original=True, zorder=zorder, clip_on=False
            )
            ax.add_collection(patch_col)  # type: ignore
    else:
        for patch in patches:
            ax.add_patch(patch)
