from __future__ import annotations

import io
import os
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import streamlit as st
from matplotlib.colors import to_hex

from pygenomeviz import __version__, load_example_dataset
from pygenomeviz.align import AlignCoord, MMseqs, MUMmer
from pygenomeviz.gui import config, plot, utils

IS_LOCAL_LAUNCH = bool(os.getenv("PGV_GUI_LOCAL"))

about_md = textwrap.dedent(
    f"""
    **pyGenomeViz v{__version__}**
    ( [GitHub](https://github.com/moshi4/pyGenomeViz) |
    [Document](https://moshi4.github.io/pyGenomeViz/) )
    """
)[1:-1]

# Streamlit page configuration
st.set_page_config(
    page_title="pyGenomeViz WebApp",
    layout="centered",
    initial_sidebar_state="expanded",
    menu_items={
        "Report a bug": "https://github.com/moshi4/pyGenomeViz/issues",
        "About": about_md,
    },
)

###########################################################
# Sidebar
###########################################################

st.sidebar.markdown(about_md)

if st.sidebar.checkbox(label="Load example genbank files", value=False):
    gbk_files = load_example_dataset("enterobacteria_phage")[0][0:4]
    gbk_list = list(map(utils.load_gbk_file, gbk_files))
else:
    with st.sidebar.expander(label="Upload User Genbank Files", expanded=True):
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
        step=5,
    )
    fig_track_height = fig_cols[1].number_input(
        label="Fig Track Height",
        value=1.0,
        min_value=0.1,
        max_value=5.0,
        step=0.5,
    )
    fig_track_ratio = fig_cols[0].number_input(
        label="Feature Track Ratio",
        value=0.3,
        min_value=0.1,
        max_value=2.0,
        step=0.1,
    )
    tick_style = fig_cols[1].selectbox(
        label="Scale",
        options=["bar", "axis", None],
        index=0,
    )
    fig_label_size = fig_cols[0].number_input(
        label="Label Size",
        value=20,
        min_value=0,
        max_value=50,
        step=5,
    )
    range_label_size = fig_cols[1].number_input(
        label="Range Label Size",
        value=15,
        min_value=0,
        max_value=50,
        step=5,
    )
    track_align_type = fig_cols[0].selectbox(
        label="Track Align Type",
        options=["left", "center", "right"],
        index=1,
    )

    fig_cfg = config.FigureConfig(
        width=fig_width,
        track_height=fig_track_height,
        tick_style=tick_style,
        align_type=str(track_align_type),
        track_ratio=fig_track_ratio,
        label_size=int(fig_label_size),
        range_label_size=int(range_label_size),
    )

with st.sidebar.expander(label="Plot Feature Options", expanded=False):
    # Feature types
    feature_types = st.multiselect(
        label="Plot Feature Types",
        options=["CDS", "Pseudo", "rRNA", "tRNA"],
        default="CDS",
        help=textwrap.dedent(
            """
            Select feature types in Genbank to be plotted  \\
            **CDS**: CDS feature with no pseudo tag  \\
            **Pseudo**: CDS feature with pseudo tag \\
            **rRNA**: rRNA feature  \\
            **tRNA**: tRNA feature
            """
        )[1:-1],
    )
    # Feature color
    color_cols = st.columns(4)
    cds_color = color_cols[0].color_picker("CDS", value=to_hex("orange"))
    pseudo_color = color_cols[1].color_picker("Pseudo", value=to_hex("grey"))
    rrna_color = color_cols[2].color_picker("rRNA", value=to_hex("lime"))
    trna_color = color_cols[3].color_picker("tRNA", value=to_hex("magenta"))
    # Feature plot style
    plotstyle_cols = st.columns(2)
    styles = ["bigarrow", "arrow", "bigbox", "box", "bigrbox", "rbox"]
    cds_plotstyle = plotstyle_cols[0].selectbox(
        label="CDS Plot Style", options=styles, index=0
    )
    pseudo_plotstyle = plotstyle_cols[1].selectbox(
        label="Pseudo Plot Style", options=styles, index=0
    )
    rrna_plotstyle = plotstyle_cols[0].selectbox(
        label="rRNA Plot Style", options=styles, index=2
    )
    trna_plotstyle = plotstyle_cols[1].selectbox(
        label="tRNA Plot Style", options=styles, index=2
    )

    label_cols = st.columns(2)
    feature_label_type = label_cols[0].selectbox(
        label="Label Type",
        options=[None, "gene", "product", "protein_id", "locus_tag"],
        index=0,
        help=textwrap.dedent(
            """
            If None, no label is shown.\\
            Otherwise, label associated with each feature in Genbank is shown.
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
        placeholder="e.g. hypothetical,unknown",
        help=textwrap.dedent(
            """
            Filter words to exclude unnecessary genbank feature labels.\\
            Multiple filter words can be set by specifying words with comma separator.\\
            (e.g. 'hypothetical,unknown')
            """
        )[1:-1],
    )
    label_filter_words = []
    for word in feature_label_filter_words.split(","):
        word = word.strip()
        if word != "":
            label_filter_words.append(word)

    show_only_top_label = st.checkbox(
        label="Show Only Top Track Label",
        value=False,
    )

    feat_cfg = config.FeatureConfig(
        types=feature_types,
        type2color=dict(
            CDS=cds_color,
            Pseudo=pseudo_color,
            rRNA=rrna_color,
            tRNA=trna_color,
        ),
        type2plotstyle=dict(
            CDS=str(cds_plotstyle),
            Pseudo=str(pseudo_plotstyle),
            rRNA=str(rrna_plotstyle),
            tRNA=str(trna_plotstyle),
        ),
        label_type=feature_label_type,
        label_size=int(feature_label_size),
        label_filter_words=label_filter_words,
        show_only_top_label=show_only_top_label,
    )

with st.sidebar.expander(label="Plot Link Options", expanded=False):
    # Check aligner installation
    aln_method_options: list[str | None] = [None]
    if MUMmer.check_installation(exit_on_false=False):
        aln_method_options.extend(["MUMmer (protein)", "MUMmer (nucleotide)"])
    if MMseqs.check_installation(exit_on_false=False):
        aln_method_options.append("MMseqs")

    aln_method = st.selectbox(
        "Genome Comparison Method",
        options=aln_method_options,
        help=textwrap.dedent(
            """
            Genome comparison method for link visualization.\\
            [MUMmer](https://github.com/mummer4/mummer) or
            [MMseqs](https://github.com/soedinglab/MMseqs2)
            installation is required to enable this functionality.

            **MUMmer (protein)**:\\
            MUMmer(promer) is used for genome alignment based on
            six frame translation between adjacent genomes\\
            **MUMmer (nucleotide)**:\\
            MUMmer(nucmer) is used for genome alignment between adjacent genomes \\
            **MMseqs**:\\
            MMseqs is used for reciprocal best-hit CDS search between adjacent genomes
            """
        )[1:-1],
    )
    link_cols = st.columns(2)
    min_length = link_cols[0].number_input(
        label="Min Length",
        value=0,
        min_value=0,
        help="Minimum length of genome comparison results to be shown",
    )
    min_identity = link_cols[1].number_input(
        label="Min Identity",
        value=0,
        min_value=0,
        max_value=100,
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
        normal_link_color=normal_link_color,
        inverted_link_color=inverted_link_color,
        curve=link_curve,
        colorbar_height=colorbar_height,
        min_length=int(min_length),
        min_identity=min_identity,
    )


###########################################################
# Main Screen
###########################################################

st.header("pyGenomeViz: Genome Visualization WebApp")

# If no genbank file exists, stop execution
if len(gbk_list) == 0:
    if not IS_LOCAL_LAUNCH:
        st.warning(
            textwrap.dedent(
                """
                :warning: This application is running on Streamlit Cloud.
                Due to the limited CPU and Memory resources of Streamlit Cloud,
                this page is intended for demonstration purposes only.
                Therefore, if you want to visualize your own genome data,
                it is recommended that you run pyGenomeViz web application
                in your local environment.
                See [pgv-gui document](https://moshi4.github.io/pyGenomeViz/gui-docs/pgv-gui/)
                for details.
                """
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

    for gbk in gbk_list:
        range_cols = st.columns([3, 3, 1])
        min_range = range_cols[0].number_input(
            label=f"**{gbk.name}** ({gbk.full_genome_length:,} bp)",
            min_value=0,
            max_value=gbk.full_genome_length,
            value=0,
            step=1,
            key=f"{gbk.name} start",
        )
        min_range = int(min_range)
        max_range = range_cols[1].number_input(
            label=f"{gbk.name} end",
            min_value=0,
            max_value=gbk.full_genome_length,
            value=gbk.full_genome_length,
            step=1,
            label_visibility="hidden",
            key=f"{gbk.name} end",
        )
        max_range = int(max_range)
        reverse = range_cols[2].selectbox(
            label="Reverse",
            options=["Yes", "No"],
            index=1,
            key=f"{gbk.name} reverse",
        )
        if min_range > max_range:
            st.error(f"**{max_range=}** must be larger than **{min_range=}**")
            st.stop()

        gbk.min_range = min_range
        gbk.max_range = max_range
        gbk.reverse = True if reverse == "Yes" else False

fig_ctl_cols = fig_ctl_container.columns([1, 2, 3])

# Plot figure
gv, fig, align_coords = plot.create_genomeviz(
    gbk_list,
    config.PgvConfig(fig_cfg, feat_cfg, aln_cfg),
)
fig_container.pyplot(fig, use_container_width=not expand_figure)

# Set figure download button
fig_format = fig_ctl_cols[0].selectbox(
    "Format",
    options=["png", "svg", "html"],
    index=0,
    format_func=str.upper,
    label_visibility="collapsed",
)

fig_save_label = f"Save Figure as {fig_format.upper()}"
filename = f"pgv_result.{fig_format}"
if fig_format in ("png", "svg"):
    fig_save_data = io.BytesIO()
    fig.savefig(fig_save_data, format=fig_format)
elif fig_format == "html":
    fig_save_data = io.BytesIO()
    gv.savefig_html(fig_save_data, fig)
else:
    raise ValueError(f"{fig_format=} is invalid.")

fig_ctl_cols[1].download_button(
    label=fig_save_label,
    data=fig_save_data,
    file_name=filename,
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
