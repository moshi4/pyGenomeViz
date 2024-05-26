from __future__ import annotations

import io
import textwrap
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import streamlit as st
from matplotlib.colors import to_hex

import pygenomeviz
from pygenomeviz.align import AlignCoord, Blast, MMseqs, MUMmer
from pygenomeviz.gui import config, plot, utils
from pygenomeviz.typing import AlnMethod
from pygenomeviz.utils import load_example_genbank_dataset

# Constant values
PLOTSTYLES = ["bigarrow", "arrow", "bigbox", "box", "bigrbox", "rbox"]
DEFAULT_FEATURE_TYPE2COLOR = defaultdict(
    lambda: to_hex("black"),
    CDS=to_hex("orange"),
    rRNA=to_hex("lime"),
    tRNA=to_hex("magenta"),
)
DEFAULT_PSEUDO_COLOR = to_hex("lightgrey")

GITHUB_URL = "https://github.com/moshi4/pyGenomeViz"
DOCS_URL = "https://moshi4.github.io/pyGenomeViz/"
ABOUD_MD = textwrap.dedent(
    f"""
    **pyGenomeViz v{pygenomeviz.__version__}**
    ( [GitHub]({GITHUB_URL}) | [Document]({DOCS_URL}) )
    """
)[1:-1]

# Streamlit page configuration
st.set_page_config(
    page_title="pyGenomeViz WebApp",
    layout="centered",
    initial_sidebar_state="expanded",
    menu_items={
        "Report a bug": "https://github.com/moshi4/pyGenomeViz/issues",
        "About": ABOUD_MD,
    },
)

###########################################################
# Sidebar
###########################################################

st.sidebar.markdown(ABOUD_MD)

checkbox_label1 = "Example Phage Genbank Files"
checkbox_label2 = "Example Bacteria Genbank Files"
if st.sidebar.checkbox(checkbox_label1):
    gbk_files = load_example_genbank_dataset("yersinia_phage")[0:4]
    gbk_list = list(map(utils.load_gbk_file, gbk_files))
elif utils.is_local_launch() and st.sidebar.checkbox(checkbox_label2):
    gbk_files = load_example_genbank_dataset("mycoplasma_mycoides")[0:4]
    gbk_list = list(map(utils.load_gbk_file, gbk_files))
else:
    with st.sidebar.expander("Upload User Genbank Files", expanded=True):
        # Genbank files upload widgets
        upload_files = st.file_uploader(
            label="Upload genbank files (\\*.gb|\\*.gbk|\\*.gbff)",
            type=["gb", "gbk", "gbff"],
            accept_multiple_files=True,
            help=(
                textwrap.dedent(
                    """
                    Genomes are displayed on each track in the order of file upload.\\
                    Genome comparison is performed between adjacent genomes.
                    """
                )[1:-1]
            ),
        )
        if upload_files is None:
            gbk_list = []
        else:
            gbk_list = list(map(utils.load_gbk_file, upload_files))

with st.sidebar.expander(label="Figure Appearance Options", expanded=False):
    fig_cols = st.columns(2)

    fig_width = fig_cols[0].number_input(
        label="Fig Width",
        value=15,
        min_value=10,
        max_value=50,
        step=1,
    )
    fig_track_height = fig_cols[1].number_input(
        label="Fig Track Height",
        value=1.0,
        min_value=0.1,
        max_value=5.0,
        step=0.1,
        format="%.1f",
    )
    feature_track_ratio = fig_cols[0].number_input(
        label="Feature Track Ratio",
        value=0.25,
        min_value=0.01,
        max_value=2.0,
        step=0.05,
    )
    link_track_ratio = fig_cols[1].number_input(
        label="Link Track Ratio",
        value=1.00,
        min_value=0.01,
        max_value=2.0,
        step=0.05,
    )
    fig_label_size = fig_cols[0].number_input(
        label="Label Size",
        value=20,
        min_value=0,
        max_value=50,
        step=1,
    )
    range_label_size = fig_cols[1].number_input(
        label="Range Label Size",
        value=0,
        min_value=0,
        max_value=20,
        step=1,
    )
    track_align_type = fig_cols[0].selectbox(
        label="Track Align Type",
        options=["left", "center", "right"],
        index=1,
    )
    scale_style = fig_cols[1].selectbox(
        label="Scale",
        options=["bar", "xticks", None],
        index=0,
    )
    seg_space_ratio = fig_cols[0].number_input(
        "Segment Space Ratio",
        value=0.02,
        min_value=0.0,
        max_value=0.1,
        step=0.01,
    )

    fig_cfg = config.FigureConfig(
        width=fig_width,
        track_height=fig_track_height,
        feature_track_ratio=feature_track_ratio,
        link_track_ratio=link_track_ratio,
        label_size=int(fig_label_size),
        range_label_size=int(range_label_size),
        track_align_type=str(track_align_type),
        scale_style=scale_style,
        seg_space_ratio=seg_space_ratio,
    )


with st.sidebar.expander(label="Plot Feature Options", expanded=False):
    feature_type2color = {}
    feature_type2plotstyle = {}
    # Feature types
    all_feature_types = utils.extract_all_feature_types(gbk_list)
    default_feature_types = []
    for feature_type in ["CDS", "rRNA"]:
        if feature_type in all_feature_types:
            default_feature_types.append(feature_type)
    select_feature_types = st.multiselect(
        label="Plot Feature Types",
        options=all_feature_types,
        default=default_feature_types,
        help=textwrap.dedent(
            """
            Genbank feature types to be shown.\\
            Users can select the feature types contained in the uploaded Genbank files.
            """
        )[1:-1],
    )
    if len(select_feature_types) >= 1:
        tab_container = st.container(border=True)
        tab_container.caption("Plotstyle & Color Options Tab")
        for feature_type, tab in zip(
            select_feature_types, tab_container.tabs(select_feature_types)
        ):
            with tab:
                cols = tab.columns([3, 1])
                plotstyle = cols[0].selectbox(
                    "Plotstyle",
                    options=PLOTSTYLES,
                    index=1,
                    key=f"{feature_type} plotstyle",
                )
                feature_type2plotstyle[feature_type] = str(plotstyle)
                color = cols[1].color_picker(
                    "Color",
                    value=DEFAULT_FEATURE_TYPE2COLOR[feature_type],
                    key=f"{feature_type} color",
                )
                feature_type2color[feature_type] = color

    cols = st.columns(2)
    line_width = cols[0].number_input(
        "Line Width",
        min_value=0.0,
        max_value=1.0,
        value=0.0,
        step=0.1,
        format="%.1f",
    )
    pseudo_color = cols[1].color_picker(
        "Pseudo Color",
        value=DEFAULT_PSEUDO_COLOR,
        help="Color of features containing '/pseudo' or '/pseudogene' tags",
    )

    label_target_track = st.radio(
        label="Label Target Track",
        options=["top", "all"],
        horizontal=True,
        format_func=lambda s: str(s).capitalize() + " Track",
    )
    label_cols = st.columns(2)
    feature_label_type = label_cols[0].selectbox(
        label="Label Type",
        options=[None, "gene", "product", "protein_id", "locus_tag"],
        index=0,
        help=textwrap.dedent(
            """
            If None, no label is shown.\\
            Otherwise, label associated with each feature is shown.
            """
        )[1:-1],
    )
    feature_label_size = label_cols[1].number_input(
        label="Label Size",
        value=10,
        min_value=5,
        max_value=20,
        step=1,
    )
    feature_label_filter_words = st.text_input(
        label="Label Filter Words",
        placeholder="e.g. unknown,putative,...",
        help=textwrap.dedent(
            """
            Filter words to exclude unnecessary genbank feature labels.\\
            Multiple filter words can be set by specifying words with comma separator (e.g. 'unknown,putative,...').\\
            Labels containing the 'hypothetical' word are excluded by default.
            """  # noqa: E501
        )[1:-1],
    )
    label_filter_words = []
    for word in feature_label_filter_words.split(","):
        word = word.strip()
        if word != "":
            label_filter_words.append(word)

    feat_cfg = config.FeatureConfig(
        types=select_feature_types,
        type2plotstyle=feature_type2plotstyle,
        type2color=feature_type2color,
        line_width=line_width,
        pseudo_color=pseudo_color,
        label_target_track=str(label_target_track),
        label_type=feature_label_type,
        label_size=int(feature_label_size),
        label_filter_words=label_filter_words,
    )

with st.sidebar.expander(label="Plot Link Options", expanded=False):
    # Check aligner installation
    aln_method_options: list[AlnMethod] = [None]
    if MUMmer.check_installation(exit_on_false=False):
        aln_method_options.extend(["MUMmer (nucleotide)", "MUMmer (protein)"])
    if MMseqs.check_installation(exit_on_false=False):
        aln_method_options.append("MMseqs RBH")
    if Blast.check_installation(exit_on_false=False):
        aln_method_options.extend(["BLAST (nucleotide)", "BLAST (protein)"])

    aln_method = st.selectbox(
        "Genome Comparison Method",
        options=aln_method_options,
        help=textwrap.dedent(
            """
            Genome comparison method between adjacent genomes for link visualization. \\
            [MUMmer](https://github.com/mummer4/mummer) or
            [MMseqs](https://github.com/soedinglab/MMseqs2) or
            [BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html)
            installation is required to enable this functionality.

            **MUMmer (nucleotide)**:\\
            MUMmer(nucmer) is used for genome alignment \\
            **MUMmer (protein)**: \\
            MUMmer(promer) is used for genome alignment based on six frame translation \\
            **MMseqs RBH**:\\
            MMseqs is used for reciprocal best-hit CDS search \\
            **BLAST (nucleotide)**: \\
            BLAST(blastn) is used for genome alignment \\
            **BLAST (protein)**: \\
            BLAST(tblastx) is used for genome alignment based on six frame translation
            """  # noqa: E501
        )[1:-1],
    )
    link_cols = st.columns(2)
    min_length = link_cols[0].number_input(
        label="Min Length",
        value=0,
        min_value=0,
        step=100,
        help="Minimum length of genome comparison results to be shown",
    )
    min_identity = link_cols[1].number_input(
        label="Min Identity",
        value=0,
        min_value=0,
        max_value=100,
        step=10,
        help="Minimum identity of genome comparison results to be shown",
    )
    link_style = link_cols[0].selectbox(
        label="Link Style",
        options=["Normal", "Curve"],
        index=0,
    )
    link_curve = True if link_style == "Curve" else False
    colorbar_height = link_cols[1].number_input(
        label="Colorbar Height",
        value=0.3,
        min_value=0.0,
        max_value=1.0,
        step=0.1,
        help=textwrap.dedent(
            """
            Colorbar is auto generated from link identity.\\
            If **0.0** is set as colorbar height, colorbar is not shown.
            """
        )[1:-1],
    )
    normal_link_color = link_cols[0].color_picker(
        label="Link (Normal)",
        value=to_hex("grey"),
    )
    inverted_link_color = link_cols[1].color_picker(
        label="Link (Inverted)",
        value=to_hex("red"),
    )

    aln_cfg = config.AlignConfig(
        method=aln_method,
        min_length=int(min_length),
        min_identity=min_identity,
        curve=link_curve,
        colorbar_height=colorbar_height,
        normal_link_color=normal_link_color,
        inverted_link_color=inverted_link_color,
    )


###########################################################
# Main Screen
###########################################################

st.header("pyGenomeViz Streamlit Web Application")

# If no genbank file exists, stop execution
if len(gbk_list) == 0:
    if not utils.is_local_launch():
        st.warning(
            textwrap.dedent(
                """
                :warning: This application is running on Streamlit Cloud.
                Due to the limited CPU and Memory constraints of Streamlit Cloud,
                this demo page limits the maximum uploadable Genbank file size to 1 MB.
                Therefore, if you want to visualize your own genome data,
                it is recommended that you run pyGenomeViz web application in your local environment.
                See [pgv-gui document](https://moshi4.github.io/pyGenomeViz/gui-docs/pgv-gui/) for details.
                """  # noqa: E501
            ),
        )
    demo_gif_file = Path(__file__).parent / "assets" / "pgv_demo.gif"
    st.image(str(demo_gif_file))
    st.stop()

expand_figure = st.checkbox(label="Expand Figure", value=False)
fig_container = st.container()
fig_ctl_container = st.container()
genome_info_container = st.container()

with genome_info_container.form(key="form"):
    title_col, form_col = st.columns([4, 1])
    title_col.markdown("**Genome Min-Max Range & Reverse Option**")
    form_col.form_submit_button(
        label="Update Figure",
        help="Apply min-max range & reverse option changes to figure",
    )

    name2seqid2range: dict[str, dict[str, tuple[int, int]]] = {}
    for gbk in gbk_list:
        expander_label = f"**{gbk.name} ({len(gbk.records)} records)**"
        with st.expander(expander_label, expanded=False):
            seqid2range = {}
            seqid2features = gbk.get_seqid2features(None)
            for idx, (seqid, size) in enumerate(gbk.get_seqid2size().items()):
                range_cols = st.columns([3, 3, 1])
                min_range = range_cols[0].number_input(
                    label=f"**{seqid}** ({size:,} bp)",
                    min_value=0,
                    max_value=size,
                    value=0,
                    step=1,
                    key=f"{gbk.name} {seqid} start",
                    help=utils.get_features_count_label(seqid2features[seqid]),
                )
                min_range = int(min_range)
                max_range = range_cols[1].number_input(
                    label=f"{gbk.name} end",
                    min_value=0,
                    max_value=size,
                    value=size,
                    step=1,
                    format="%d",
                    label_visibility="hidden",
                    key=f"{gbk.name} {seqid} end",
                )
                max_range = int(max_range)
                if min_range > max_range:
                    st.error(f"**{max_range=}** must be larger than **{min_range=}**")
                    st.stop()
                if min_range != max_range:
                    seqid2range[seqid] = (min_range, max_range)
                reverse = range_cols[2].selectbox(
                    label="Reverse",
                    options=[True, False],
                    index=1,
                    format_func=lambda b: "Yes" if b else "No",
                    key=f"{gbk.name} {seqid} reverse",
                )
                if reverse is True:
                    gbk.records[idx] = gbk.records[idx].reverse_complement(
                        id=True, name=True, description=True
                    )
            name2seqid2range[gbk.name] = seqid2range


fig_ctl_cols = fig_ctl_container.columns([1, 2, 3])

# Plot figure
gv, align_coords = plot.plot_by_gui_cfg(
    gbk_list,
    config.PgvGuiPlotConfig(fig_cfg, feat_cfg, aln_cfg, name2seqid2range),
)
fig = gv.plotfig()
fig_container.pyplot(fig, use_container_width=not expand_figure)

# Set figure download button
fig_format = fig_ctl_cols[0].selectbox(
    "Format",
    options=["png", "svg", "html"],
    index=0,
    format_func=str.upper,
    label_visibility="collapsed",
)

fig_bytes = io.BytesIO()
if fig_format in ("png", "svg"):
    fig.savefig(fig_bytes, format=fig_format)
elif fig_format == "html":
    gv.savefig_html(fig_bytes)
else:
    raise ValueError(f"{fig_format=} is invalid.")

fig_ctl_cols[1].download_button(
    label=f"Save Figure as {fig_format.upper()}",
    data=fig_bytes,
    file_name=f"pgv_result.{fig_format}",
)

if align_coords:
    comparison_result_data = io.BytesIO()
    AlignCoord.write(align_coords, comparison_result_data)
    fig_ctl_cols[2].download_button(
        label="Save Comparison Result",
        data=comparison_result_data,
        file_name="pgv_comparison_result.tsv",
    )

# Clear & close figure to suppress memory leak
fig.clear()
plt.close(fig)
