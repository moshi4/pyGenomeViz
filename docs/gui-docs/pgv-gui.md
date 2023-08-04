# pgv-gui GUI Document

`pgv-gui` is a command to launch the GUI (Web browser) version of pyGenomeViz.
Users can easily perform data visualization of Genbank files and genome comparison
results using MUMmer or MMseqs.

![pygenomeviz_gui.gif](https://raw.githubusercontent.com/moshi4/pyGenomeViz/main/src/pygenomeviz/gui/assets/pgv_demo.gif)

## Installation

Additional installation of [streamlit](https://github.com/streamlit/streamlit) is required.
MUMmer and MMseqs are also required to enable the genome comparison functionality.

### Conda

    conda install -c conda-forge -c bioconda pygenomeviz streamlit mummer mmseqs2

### Pip

    pip install pygenomeviz[all]

In Ubuntu22.04, MUMmer and MMseqs can be installed with apt command below.

    sudo apt install mummer mmseqs

### Docker

    docker run --rm -p 8501:8501 ghcr.io/moshi4/pygenomeviz:latest pgv-gui

## Usage

The following command launches a GUI browser, which can be accessed at <https://localhost:8501>.

    pgv-gui
