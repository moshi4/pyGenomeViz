# pgv-gui Web Application Document

`pgv-gui` command is used to launch the pyGenomeViz web application.
It is developed with the streamlit web application framework,
and users can easily visualize the genome data of Genbank files and their comparison results with GUI.

![pygenomeviz_gui.gif](https://raw.githubusercontent.com/moshi4/pyGenomeViz/main/src/pygenomeviz/gui/assets/pgv_demo.gif)
**Fig.pyGenomeViz web application example ([Demo Page](https://pygenomeviz.streamlit.app))**

## Installation

Additional installation of [streamlit](https://github.com/streamlit/streamlit) is required.
MUMmer and MMseqs are also required to enable the genome comparison functionality.

### Conda

    conda install -c conda-forge -c bioconda pygenomeviz streamlit mummer mmseqs2

### Pip

    pip install pygenomeviz[gui]

In Ubuntu22.04, MUMmer and MMseqs can be installed with apt command below.

    sudo apt install mummer mmseqs2

### Docker

    docker run --rm -p 8501:8501 ghcr.io/moshi4/pygenomeviz:latest pgv-gui

## Usage

The following command launches web application, which can be accessed at <http://localhost:8501>.

    pgv-gui

By uploading the user's Genbank files in the browser, a visualization figure of the genome data is automatically displayed.
By changing the value of each widget on the browser, user can adjust the appearance of the figure,
change the genome comparison method, etc. interactively.
