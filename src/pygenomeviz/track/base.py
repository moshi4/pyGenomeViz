from __future__ import annotations

from matplotlib.axes import Axes


class Track:
    """Track Base Class"""

    def __init__(
        self,
        name: str,
        *,
        ratio: float = 1.0,
        zorder: float = 0.0,
    ):
        self._name = name
        self._ratio = ratio
        self._zorder = zorder
        self._xlim: tuple[int, int] | None = None
        self._ylim: tuple[float, float] = (-1.0, 1.0)
        self._ax: Axes | None = None

    ############################################################
    # Property
    ############################################################

    @property
    def name(self) -> str:
        """Track name"""
        return self._name

    @property
    def ratio(self) -> float:
        """Track size ratio"""
        return self._ratio

    @property
    def zorder(self) -> float:
        """Track zorder"""
        return self._zorder

    @property
    def xlim(self) -> tuple[int, int]:
        """Track axes x min-max tuple"""
        if self._xlim is None:
            raise ValueError("'xlim' is not defined!!")
        return self._xlim

    @property
    def ylim(self) -> tuple[float, float]:
        """Track axes y min-max tuple"""
        return self._ylim

    @property
    def ax(self) -> Axes:
        """Track axes

        Can't access ax property before calling GenomeViz class `plotfig` method.

        Returns
        -------
        ax : Axes
            Matplotlib axes
        """
        if self._ax is None:
            err_msg = "Can't access ax property before calling 'plotfig' method."
            raise ValueError(err_msg)
        return self._ax

    ############################################################
    # Public Method
    ############################################################

    def set_ratio(self, ratio: float) -> None:
        """Set track size ratio"""
        if ratio < 0:
            raise ValueError(f"{ratio=} is invalid (Must be 'ratio >= 0').")
        self._ratio = ratio

    def set_xlim(self, xlim: tuple[int, int]) -> None:
        """Set track xlim"""
        self._xlim = xlim

    def set_ax(self, ax: Axes, show_axis: bool = False) -> None:
        """Set track axes

        Parameters
        ----------
        ax : Axes
            Matplotlib axes
        show_axis : bool, optional
            Show axis for debug purpose
        """
        self._initalize_axes(ax, show_axis)
        self._ax = ax

    ############################################################
    # Private Method
    ############################################################

    def _initalize_axes(self, ax: Axes, show_axis: bool = False) -> None:
        """Initialize axes properties

        Parameters
        ----------
        ax : Axes
            Matplotlib axes
        show_axis : bool, optional
            Show axis for debug purpose
        """
        # Set spines visibility
        for pos in ("left", "right", "top", "bottom"):
            ax.spines[pos].set_visible(show_axis)
        # Set ticks visibility
        ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
        # Set zorder
        ax.set_zorder(self.zorder)
        # Set xlim
        ax.set_xlim(self.xlim)
        ax.set_ylim(*self.ylim)
        # Set facecolor to transparent
        ax.set_facecolor("none")
