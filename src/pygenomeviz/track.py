from __future__ import annotations

import math
from typing import Any, Callable, Literal

from Bio.SeqFeature import SeqFeature
from matplotlib.axes import Axes

from pygenomeviz.config import LiteralTypes
from pygenomeviz.feature import ExonFeature, Feature
from pygenomeviz.link import Link
from pygenomeviz.parser import Genbank, Gff


class Track:
    """Track BaseClass"""

    def __init__(
        self,
        name: str,
        size: int,
        labelsize: int = 0,
        spines: bool = False,
        ratio: float = 1.0,
    ):
        """
        Parameters
        ----------
        name : str
            Track name
        size : int
            Track size
        labelsize : int, optional
            Track label size
        spines : bool, optional
            Display track spines
        ratio : float, optional
            Track height ratio
        """
        self.name = name
        self.size = size
        self.labelsize = labelsize
        self.spines = spines
        self.ratio = ratio
        self._offset: int | None = None
        self._ax: Axes | None = None

    @property
    def zorder(self) -> float:
        """Track zorder"""
        return 0

    @property
    def ylim(self) -> tuple[float, float]:
        """Track y min-max limit tuple"""
        return (-1.0, 1.0)

    @property
    def tick_params(self) -> dict[str, bool]:
        """Track tick parameters dict"""
        return {
            "left": False,
            "labelleft": False,
            "bottom": False,
            "labelbottom": False,
        }

    @property
    def spines_params(self) -> dict[str, bool]:
        """Spines parameters dict"""
        return {
            "left": self.spines,
            "right": self.spines,
            "top": self.spines,
            "bottom": self.spines,
        }

    @property
    def offset(self) -> int:
        """Track alignment offset"""
        if self._offset is None:
            err_msg = "Can't access offset property before calling 'plotfig' method."
            raise ValueError(err_msg)
        return self._offset

    @property
    def ax(self) -> Axes:
        """Track axes

        Can't access ax property before calling GenomeViz class `plotfig` method.

        Returns
        -------
        ax : Axes
            Track matplotlib axes object
            (Default: xlim=`(0, gv.max_track_size)`, ylim=`(-1, 1)`)
        """
        if self._ax is None:
            err_msg = "Can't access ax property before calling 'plotfig' method."
            raise ValueError(err_msg)
        return self._ax


class FeatureTrack(Track):
    """FeatureTrack Class"""

    def __init__(
        self,
        name: str,
        size: int,
        start_pos: int = 0,
        labelsize: int = 20,
        labelcolor: str = "black",
        labelmargin: float = 0.01,
        linewidth: int = 1,
        linecolor: str = "black",
        spines: bool = False,
        ratio: float = 1.0,
    ):
        """
        Parameters
        ----------
        name : str
            Track name
        size : int
            Track size
        start_pos : int, optional
            Track start position.
            Track start-end range is defined as (start_pos, start_pos + size).
        labelsize : int, optional
            Track label size
        labelcolor : str, optional
            Track label color
        labelmargin : float, optional
            Track label margin
        linewidth : int, optional
            Track line width
        linecolor : str, optional
            Track line color
        spines : bool, optional
            Display track spines
        ratio : float, optional
            Track height ratio
        """
        super().__init__(name, size, labelsize, spines, ratio)
        self._start_pos = start_pos
        self.labelcolor = labelcolor
        self.labelmargin = labelmargin
        self.linewidth = linewidth
        self.linecolor = linecolor
        self.features: list[Feature] = []
        self.subtracks: list[FeatureSubTrack] = []

        # Sublabel parameters
        self._sublabel_y: float | None = None
        self._sublabel_text: str | None = None
        self._sublabel_size: int = 10
        self._sublabel_color: str = "black"
        self._sublabel_ha: str = "left"
        self._sublabel_va: str = "top"
        self._sublabel_kws: dict[str, Any] = {}

    @property
    def start(self) -> int:
        """Track start position"""
        return self._start_pos

    @property
    def end(self) -> int:
        """Track end position"""
        return self.start + self.size

    @property
    def zorder(self) -> float:
        """Track zorder"""
        return 10

    @property
    def label_params(self) -> dict[str, Any]:
        """Label drawing parameters"""
        return {
            "s": self.name,
            "fontsize": self.labelsize,
            "color": self.labelcolor,
            "ha": "right",
            "va": "center",
            "gid": f"TrackLabel_{self.name}",
        }

    @property
    def sublabel_params(self) -> dict[str, Any]:
        """Sublabel drawing parameters"""
        return {
            "y": self._sublabel_y,
            "s": self._sublabel_text,
            "fontsize": self._sublabel_size,
            "color": self._sublabel_color,
            "ha": self._sublabel_ha,
            "va": self._sublabel_va,
            "zorder": 10,
            **self._sublabel_kws,
        }

    def _within_valid_range(
        self, start: int, end: int, raise_error_on_false: bool = False
    ) -> bool:
        """Check if start-end position within valid track range

        Parameters
        ----------
        start : int
            Start position
        end : int
            End position
        raise_error_on_false : bool, optional
            If start-end range is invalid in track, raise error

        Returns
        -------
        check_result : bool
            Check result
        """
        if self.start <= start <= self.end and self.start <= end <= self.end:
            return True
        else:
            if raise_error_on_false:
                err_msg = f"'{self.name}' track start-end range must be "
                err_msg += f"'{self.start} <= start <= end <= {self.end}' "
                err_msg += f"({start=}, {end=} is invalid)"
                raise ValueError(err_msg)
            return False

    def add_subtrack(
        self, name: str | None = None, ratio: float = 1.0, position: str = "below"
    ) -> None:
        """Add subtrack to feature track

        Parameters
        ----------
        name : str | None, optional
            Subtrack name. If None, subtrack name is automatically set.
        ratio : float, optional
            Subtrack size ratio to feature track
        position : str, optional
            Subtrack position (`above`|`below`)
        """
        subtrack_ratio = self.ratio * ratio
        subtrack_idx = len(self.subtracks) + 1
        if name is None:
            name = f"{self.name}_subtrack{subtrack_idx}"

        if name in [t.name for t in self.subtracks]:
            raise ValueError(f"subtrack.name='{name}' is already exists.")
        if position not in ("above", "below"):
            raise ValueError(f"{position=} is invalid ('above'|'below').")

        subtrack = FeatureSubTrack(
            name, self.size, self.spines, subtrack_ratio, position
        )
        self.subtracks.append(subtrack)

    def get_subtrack(self, name: str) -> FeatureSubTrack:
        """Get subtrack by name

        Parameters
        ----------
        name : str
            Target subtrack name

        Returns
        -------
        subtrack : FeatureSubTrack
            Target subtrack
        """
        name2subtrack = {t.name: t for t in self.subtracks}
        if name not in name2subtrack.keys():
            err_msg = f"subtrack.name='{name}' is not found."
            raise ValueError(err_msg)
        return name2subtrack[name]

    def set_sublabel(
        self,
        text: str | None = None,
        size: int = 15,
        color: str = "black",
        position: str = "bottom-left",
        ymargin: float = 0.2,
        sublabel_kws: dict[str, Any] | None = None,
    ) -> None:
        """Set sublabel to feature track

        Parameters
        ----------
        text : str | None, optional
            Sublabel text. If None, `{start} - {end} bp` label text is set.
        size : int, optional
            Sublabel size
        color : str, optional
            Sublabel color
        position : str, optional
            Sublabel position (`[top|bottom]-[left|center|right]`)
        ymargin : float, optional
            Sublabel y-margin
        sublabel_kws: dict[str, Any] | None, optional
            Text properties (e.g. `dict(size=12, color="red", ...)`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.text.html>
        """
        # Sublabel text setting (e.g. '0 - 1000 bp')
        default_text = f"{self.start} - {self.end} bp"
        self._sublabel_text = default_text if text is None else text
        self._sublabel_size = size
        self._sublabel_color = color
        self._sublabel_kws = {} if sublabel_kws is None else sublabel_kws
        # Sublabel position setting
        vpos, hpos = position.split("-")
        vpos_types, hpos_types = ("top", "bottom"), ("left", "center", "right")
        if vpos not in vpos_types or hpos not in hpos_types:
            err_msg = f"{position=} is invalid pattern. "
            err_msg += "position must be '[top|bottom]-[left|center|right]'"
            raise ValueError(err_msg)
        vpos2y = {"top": 1 + ymargin, "bottom": -1 - ymargin}
        vpos2va = {"top": "bottom", "bottom": "top"}
        self._sublabel_y = vpos2y[vpos]
        self._sublabel_va = vpos2va[vpos]
        self._sublabel_ha = hpos

    def add_feature(
        self,
        start: int,
        end: int,
        strand: LiteralTypes.STRAND = 1,
        label: str = "",
        labelsize: int = 15,
        labelcolor: str = "black",
        plotstyle: LiteralTypes.PLOTSTYLE = "bigarrow",
        facecolor: str = "orange",
        edgecolor: str = "black",
        linewidth: float = 0,
        labelrotation: int = 45,
        labelvpos: LiteralTypes.LABELVPOS = "strand",
        labelhpos: LiteralTypes.LABELHPOS = "center",
        labelha: LiteralTypes.LABELHA = "left",
        arrow_shaft_ratio: float = 0.5,
        size_ratio: float = 1.0,
        patch_kws: dict[str, Any] | None = None,
    ) -> None:
        """Add feature to track

        Parameters
        ----------
        start : int
            Featrue start position
        end : int
            Feature end position
        strand : int, optional
            Feature strand
        label : str, optional
            Feature label
        labelsize : int, optional
            Feature label size
        labelcolor : str, optional
            Feature label color
        plotstyle : str, optional
            Feature plot style (`bigarrow`|`arrow`|`bigbox`|`box`|`bigrbox`|`rbox`)
        facecolor : str, optional
            Feature facecolor
        edgecolor : str, optional
            Feature edgecolor
        linewidth : float, optional
            Feature edge linewidth
        labelrotation : int, optional
            Feature label rotation
        labelvpos : str, optional
            Feature label vertical position (`top`|`center`|`bottom`|`strand`)
            If 'strand' is set, 'top' or 'bottom' is auto selected by strand.
        labelhpos : str, optional
            Feature label horizontal position (`left`|`center`|`right`)
        labelha : str, optional
            Feature label horizontal alignment (`left`|`center`|`right`)
        arrow_shaft_ratio : float, optional
            Feature arrow shaft ratio
        size_ratio : float, optional
            Feature size ratio to track
        patch_kws : dict[str, Any] | None, optional
            Patch properties (e.g. `fc="red", ec="black", lw=0.5, ...`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        # Check if start & end positions are within appropriate track range
        self._within_valid_range(start, end, raise_error_on_false=True)

        self.features.append(
            Feature(
                start,
                end,
                strand,
                label,
                labelsize,
                labelcolor,
                plotstyle,
                facecolor,
                edgecolor,
                linewidth,
                labelrotation,
                labelvpos,
                labelhpos,
                labelha,
                arrow_shaft_ratio,
                size_ratio,
                patch_kws,
            )
        )

    def add_exon_feature(
        self,
        exon_regions: list[tuple[int, int]],
        strand: LiteralTypes.STRAND = 1,
        label: str = "",
        labelsize: int = 15,
        labelcolor: str = "black",
        plotstyle: LiteralTypes.PLOTSTYLE = "bigarrow",
        facecolor: str = "orange",
        edgecolor: str = "black",
        linewidth: float = 0,
        labelrotation: int = 45,
        labelvpos: LiteralTypes.LABELVPOS = "strand",
        labelhpos: LiteralTypes.LABELHPOS = "center",
        labelha: LiteralTypes.LABELHA = "left",
        arrow_shaft_ratio: float = 0.5,
        size_ratio: float = 1.0,
        exon_labels: list[str] | None = None,
        exon_label_kws: dict[str, Any] | None = None,
        patch_kws: dict[str, Any] | None = None,
        intron_patch_kws: dict[str, Any] | None = None,
    ) -> None:
        """Add exon feature to track

        Parameters
        ----------
        exon_regions : list[tuple[int, int]]
            Exon feature start-end postion list
        strand : int, optional
            Feature strand
        label : str, optional
            Feature label
        labelsize : int, optional
            Feature label size
        labelcolor : str, optional
            Feature label color
        plotstyle : str, optional
            Feature plot style (`bigarrow`|`arrow`|`bigbox`|`box`|`bigrbox`|`rbox`)
        facecolor : str, optional
            Feature facecolor
        edgecolor : str, optional
            Feature edgecolor
        linewidth : float, optional
            Feature edge linewidth
        labelrotation : int, optional
            Feature label rotation
        labelvpos : str, optional
            Feature label vertical position (`top`|`center`|`bottom`|`strand`)
            If 'strand' is set, 'top' or 'bottom' is auto selected by strand.
        labelhpos : str, optional
            Feature label horizontal position (`left`|`center`|`right`)
        labelha : str, optional
            Feature label horizontal alignment (`left`|`center`|`right`)
        arrow_shaft_ratio : float, optional
            Feature arrow shaft ratio
        size_ratio : float, optional
            Feature size ratio to track
        exon_labels: list[str] | None, optional
            Exon labels. Array length must be same as `exon_regions`.
        exon_label_kws : dict[str, Any] | None, optional
            Text properties (e.g. `dict(size=12, color="red", ...)`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.text.html>
        patch_kws : dict[str, Any] | None, optional
            Exon patch properties (e.g. `dict(fc="red", ec="black", lw=0.5, ...)`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        intron_patch_kws : dict[str, Any] | None, optional
            Intron patch properties (e.g. `dict(ec="red", lw=2, ...)`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        # Check if start & end positions are within appropriate track range
        for (exon_start, exon_end) in exon_regions:
            self._within_valid_range(exon_start, exon_end, raise_error_on_false=True)

        self.features.append(
            ExonFeature(
                exon_regions,
                strand,
                label,
                labelsize,
                labelcolor,
                plotstyle,
                facecolor,
                edgecolor,
                linewidth,
                labelrotation,
                labelvpos,
                labelhpos,
                labelha,
                arrow_shaft_ratio,
                size_ratio,
                exon_labels,
                exon_label_kws,
                patch_kws,
                intron_patch_kws,
            )
        )

    def add_genbank_features(
        self,
        gbk: Genbank,
        feature_type: str = "CDS",
        label_type: str | None = None,
        label_handle_func: Callable[[str], str] | None = None,
        allow_partial: bool = False,
        pseudogene: bool = False,
        labelsize: int = 15,
        labelcolor: str = "black",
        plotstyle: LiteralTypes.PLOTSTYLE = "bigarrow",
        facecolor: str = "orange",
        facecolor_handle_func: Callable[[SeqFeature], str] | None = None,
        edgecolor: str = "black",
        linewidth: float = 0,
        labelrotation: int = 45,
        labelvpos: LiteralTypes.LABELVPOS = "strand",
        labelhpos: LiteralTypes.LABELHPOS = "center",
        labelha: LiteralTypes.LABELHA = "left",
        arrow_shaft_ratio: float = 0.5,
        size_ratio: float = 1.0,
        patch_kws: dict[str, Any] | None = None,
    ) -> None:
        """Add features from genbank record

        Parameters
        ----------
        gbk : Genbank
            Genbank object
        feature_type : str, optional
            Feature type (e.g. `CDS`,`rRNA`,`tRNA`,etc...)
        label_type : str | None, optional
            Label type (e.g. `gene`,`protein_id`,`product`,etc...)
        label_handle_func : Callable[[str], str] | None, optional
            User defined function to handle label.
            Useful for filtering out unnecesary labels such as `hypothetical ~~~`,
            omitting labels with long characters, etc.
        allow_partial : bool, optional
            If True, features that are partially included in range are also extracted
        pseudogene : bool, optional
            If True and `feature_type='CDS'`, only add CDS features with
            `/pseudo` or `/pseudogene` qualifiers.
        labelsize : int, optional
            Feature label size
        labelcolor : str, optional
            Feature label color
        plotstyle : str, optional
            Feature plot style (`bigarrow`|`arrow`|`bigbox`|`box`|`bigrbox`|`rbox`)
        facecolor : str, optional
            Feature facecolor.
            If Genbank qualifiers has facecolor key (e.g. `/facecolor="red"`),
            facecolor key value is applied preferentially.
        facecolor_handle_func : Callable[[SeqFeature], str] | None, optional
            User defined function to handle feature facecolor.
        edgecolor : str, optional
            Feature edgecolor
        linewidth : float, optional
            Feature edge linewidth
        labelrotation : int, optional
            Feature label rotation
        labelvpos : str, optional
            Feature label vertical position (`top`|`center`|`bottom`|`strand`)
            If 'strand' is set, 'top' or 'bottom' is auto selected by strand.
        labelhpos : str, optional
            Feature label horizontal position (`left`|`center`|`right`)
        labelha : str, optional
            Feature label horizontal alignment (`left`|`center`|`right`)
        arrow_shaft_ratio : float, optional
            Feature arrow shaft ratio
        size_ratio : float, optional
            Feature size ratio to track
        patch_kws : dict[str, Any] | None, optional
            Patch properties (e.g. `dict(fc="red", ec="black", lw=0.5, ...)`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        target_features = gbk.extract_features(
            feature_type, None, False, allow_partial, pseudogene
        )
        for feature in target_features:
            # Check if start & end positions are within appropriate track range
            start = int(str(feature.location.start))
            end = int(str(feature.location.end))
            self._within_valid_range(start, end, raise_error_on_false=True)

            # Get label value & apply handle func if exists
            label = feature.qualifiers.get(label_type, [""])[0]
            if label_handle_func is not None:
                label = label_handle_func(label)
            # Get facecolor & apply handle func if exists
            facecolor = feature.qualifiers.get("facecolor", [facecolor])[0]
            if facecolor_handle_func is not None:
                facecolor = facecolor_handle_func(feature)

            self.features.append(
                Feature(
                    start,
                    end,
                    feature.strand,
                    label,
                    labelsize,
                    labelcolor,
                    plotstyle,
                    facecolor,
                    edgecolor,
                    linewidth,
                    labelrotation,
                    labelvpos,
                    labelhpos,
                    labelha,
                    arrow_shaft_ratio,
                    size_ratio,
                    patch_kws,
                    feature,
                )
            )

    def add_gff_features(
        self,
        gff: Gff,
        feature_type: str = "CDS",
        parse_exon_intron: bool = False,
        label_type: str | None = None,
        label_handle_func: Callable[[str], str] | None = None,
        labelsize: int = 15,
        labelcolor: str = "black",
        plotstyle: LiteralTypes.PLOTSTYLE = "bigarrow",
        facecolor: str = "orange",
        facecolor_handle_func: Callable[[SeqFeature], str] | None = None,
        edgecolor: str = "black",
        linewidth: float = 0,
        labelrotation: int = 45,
        labelvpos: LiteralTypes.LABELVPOS = "strand",
        labelhpos: LiteralTypes.LABELHPOS = "center",
        labelha: LiteralTypes.LABELHA = "left",
        arrow_shaft_ratio: float = 0.5,
        size_ratio: float = 1.0,
        patch_kws: dict[str, Any] | None = None,
        intron_patch_kws: dict[str, Any] | None = None,
    ) -> None:
        """Add features from GFF record

        Parameters
        ----------
        gff : Gff
            GFF object
        feature_type : str, optional
            Feature type (e.g. `CDS`,`gene`,`mRNA`,etc...)
        parse_exon_intron : bool, optional
            If True, try to parse and add exon-intron structured feature
            (Expected to be used for eukaryote with feature_type=`mRNA`|`ncRNA`)
        label_type : str | None, optional
            Label type in attributes column (e.g. `ID`,`Name`,`product`,etc...)
        label_handle_func : Callable[[str], str] | None, optional
            User defined function to handle label.
            Useful for filtering out unnecesary labels such as `hypothetical ~~~`,
            omitting labels with long characters, etc.
        labelsize : int, optional
            Feature label size
        labelcolor : str, optional
            Feature label color
        plotstyle : str, optional
            Feature plot style (`bigarrow`|`arrow`|`bigbox`|`box`|`bigrbox`|`rbox`)
        facecolor : str, optional
            Feature facecolor.
            If GFF attributes column has facecolor tag (e.g. `facecolor=red;`),
            facecolor tag value is applied preferentially.
        facecolor_handle_func : Callable[[SeqFeature], str] | None, optional
            User defined function to handle feature facecolor.
        edgecolor : str, optional
            Feature edgecolor
        linewidth : float, optional
            Feature edge linewidth
        labelrotation : int, optional
            Feature label rotation
        labelvpos : str, optional
            Feature label vertical position (`top`|`center`|`bottom`|`strand`)
            If 'strand' is set, 'top' or 'bottom' is auto selected by strand.
        labelhpos : str, optional
            Feature label horizontal position (`left`|`center`|`right`)
        labelha : str, optional
            Feature label horizontal alignment (`left`|`center`|`right`)
        arrow_shaft_ratio : float, optional
            Feature arrow shaft ratio
        size_ratio : float, optional
            Feature size ratio to track
        patch_kws : dict[str, Any] | None, optional
            Patch properties (e.g. `dict(fc="red", ec="black", lw=0.5, ...)`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        intron_patch_kws : dict[str, Any] | None, optional
            Intron patch properties (e.g. `dict(ec="red", lw=0.5, ...)`)
            <https://matplotlib.org/stable/api/_as_gen/matplotlib.patches.Patch.html>
        """
        # Extract target features
        if parse_exon_intron:
            features = gff.extract_exon_features(feature_type)
        else:
            features = gff.extract_features(feature_type)

        # Add features in track
        for feature in features:
            # Get label value & apply handle func if exists
            label = feature.qualifiers.get(label_type, [""])[0]
            if label_handle_func is not None:
                label = label_handle_func(label)
            # Get facecolor & apply handle func if exists
            facecolor = feature.qualifiers.get("facecolor", [facecolor])[0]
            if facecolor_handle_func is not None:
                facecolor = facecolor_handle_func(feature)
            # Create feature plot property dict
            feature_kws = dict(
                strand=feature.strand,
                label=label,
                labelsize=labelsize,
                labelcolor=labelcolor,
                plotstyle=plotstyle,
                facecolor=facecolor,
                edgecolor=edgecolor,
                linewidth=linewidth,
                labelrotation=labelrotation,
                labelvpos=labelvpos,
                labelhpos=labelhpos,
                labelha=labelha,
                arrow_shaft_ratio=arrow_shaft_ratio,
                size_ratio=size_ratio,
                patch_kws=patch_kws,
                seq_feature=feature,
            )

            if parse_exon_intron:
                # Extract exon regions
                exon_regions: list[tuple[int, int]] = []
                for loc in feature.location.parts:
                    start, end = int(str(loc.start)), int(str(loc.end))
                    self._within_valid_range(start, end, raise_error_on_false=True)
                    exon_regions.append((start, end))
                # Add exon feature
                feature_kws["intron_patch_kws"] = intron_patch_kws
                self.features.append(ExonFeature(exon_regions, **feature_kws))
            else:
                start = int(str(feature.location.start))
                end = int(str(feature.location.end))
                self._within_valid_range(start, end, raise_error_on_false=True)
                # Add feature
                self.features.append(Feature(start, end, **feature_kws))

    def __str__(self):
        return f"{self.name}:{self.start}-{self.end}"


class FeatureSubTrack(Track):
    """FeatureSubTrack Class"""

    def __init__(
        self,
        name: str,
        size: int,
        spines: bool = False,
        ratio: float = 1.0,
        position: Literal["above", "below"] = "below",
    ):
        """
        Parameters
        ----------
        name : str
            Subtrack name
        size : int
            Subtrack size
        spines : bool, optional
            Display track spines
        ratio : float, optional
            Subtrack height ratio
        position : str, optional
            Subtrack position (`above`|`below`)
        """
        super().__init__(name, size, 0, spines, ratio)
        self.position = position


class LinkTrack(Track):
    """LinkTrack Class"""

    def __init__(self, name: str, spines: bool = False, ratio: float = 1.0):
        """
        Parameters
        ----------
        name : str
            Track name
        spines : bool, optional
            Display track spines
        ratio : float, optional
            Track height ratio
        """
        super().__init__(name, 0, 0, spines, ratio)
        self.links: list[Link] = []

    def add_link(self, link: Link) -> None:
        """Add link to track

        Parameters
        ----------
        link : Link
            Link between feature track
        """
        self.links.append(link)


class TickTrack(Track):
    """TickTrack Class"""

    def __init__(
        self,
        size: int,
        labelsize: int = 15,
        spines: bool = False,
        ratio: float = 1.0,
        tick_style: str = "bar",
    ):
        """
        Parameters
        ----------
        size : int
            Track size
        labelsize : int, optional
            Tick label size
        spines : bool, optional
            Display track spines
        ratio : float, optional
            Track height ratio
        tick_style : str, optional
            Tick style (`axis`|`bar`)
        """
        super().__init__("tick", size, labelsize, spines, ratio)
        self.tick_style = tick_style
        self.height_scale = 1.0

    @property
    def tick_params(self) -> dict[str, Any]:
        """Track tick parameters dict"""
        return {
            "left": False,
            "labelleft": False,
            "bottom": self.tick_style == "axis",
            "labelbottom": self.tick_style == "axis",
            "labelsize": self.labelsize,
        }

    @property
    def spines_params(self) -> dict[str, bool]:
        """Spines parameters dict"""
        return {
            "left": self.spines,
            "right": self.spines,
            "top": self.spines,
            "bottom": self.spines or self.tick_style == "axis",
        }

    @property
    def scalebar_text_params(self) -> dict[str, Any]:
        """Scalebar text parameters dict"""
        return {
            "x": self.xcenter,
            "y": self.ymin,
            "s": self.scalebar_label,
            "fontsize": self.labelsize,
            "ha": "center",
            "va": "top",
        }

    @property
    def unit(self) -> str:
        """Unit (bp, Kb, Mb, Gb)"""
        for unit, value in self.unit2base_value.items():
            if self.size >= value:
                return unit
        raise ValueError("Unexpected error.")

    @property
    def format_str(self) -> str:
        """Format string ('.0f' or '.1f')"""
        if self.size / self.base_value >= 10:
            return ".0f"
        else:
            return ".1f"

    @property
    def base_value(self) -> float:
        """Base value"""
        return self.unit2base_value[self.unit]

    @property
    def unit2base_value(self) -> dict[str, int]:
        """Unit & base value dict"""
        return {"Gb": 10**9, "Mb": 10**6, "Kb": 10**3, "bp": 1}

    @property
    def scalebar_label(self) -> str:
        """Label"""
        label = f"{self.scalebar_size / self.base_value:{self.format_str}}"
        return f"{label} {self.unit}"

    @property
    def xmin(self) -> float:
        """xmin"""
        return self.size - self.scalebar_size

    @property
    def xmax(self) -> float:
        """xmax"""
        return self.size

    @property
    def xcenter(self) -> float:
        """xcenter"""
        return (self.xmin + self.xmax) / 2

    @property
    def ymin(self) -> float:
        """ymin"""
        return self.ylim[0] - abs(self.ylim[0] * 0.1 * self.height_scale)

    @property
    def ycenter(self) -> float:
        """ycenter"""
        return self.ylim[0]

    @property
    def ymax(self) -> float:
        """ymax"""
        return self.ylim[0] + abs(self.ylim[0] * 0.1 * self.height_scale)

    @property
    def scalebar_size(self) -> float:
        """Scalebar size"""
        min_scalebar_size = self.size * 0.1
        unit = int(10 ** (len(str(int(min_scalebar_size))) - 1))
        value = math.ceil(min_scalebar_size / unit)
        steps = [1, 2, 5, 10]
        for i in range(0, len(steps) - 1):
            if steps[i] < value < steps[i + 1]:
                return steps[i + 1] * unit
        return value * unit

    def tick_formatter(self, value: float, pos: int) -> str:
        """Tick formatter

        Parameters
        ----------
        value : float
            Format target tick value
        pos : int
            Tick position (Not used for value formatting)

        Returns
        -------
        format_value : str
            Format tick value string
        """
        tick_value = value / self.base_value
        return f"{tick_value:{self.format_str}} {self.unit}"
