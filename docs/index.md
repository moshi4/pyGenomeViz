# pyGenomeViz

![Python3](https://img.shields.io/badge/Language-Python3-steelblue)
![OS](https://img.shields.io/badge/OS-_Windows_|_Mac_|_Linux-steelblue)
![License](https://img.shields.io/badge/License-MIT-steelblue)
[![Latest PyPI version](https://img.shields.io/pypi/v/pygenomeviz.svg)](https://pypi.python.org/pypi/pygenomeviz)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/pygenomeviz.svg?color=green)](https://anaconda.org/bioconda/pygenomeviz)  

## Overview

pyGenomeViz is a genome visualization python package for comparative genomics implemented based on matplotlib.
This package is developed for the purpose of easily and beautifully plotting genomic
features and sequence similarity comparison links between multiple genomes.
It supports genome visualization of Genbank/GFF format file and can be saved figure in various formats (JPG/PNG/SVG/PDF/HTML).
User can use pyGenomeViz for interactive genome visualization figure plotting on jupyter notebook,
or automatic genome visualization figure plotting in genome analysis scripts/pipelines.

<figure markdown>
  ![pygenomeviz_gallery.png](./images/pygenomeviz_gallery.png)
  <figcaption>pyGenomeViz example plot gallery</figcaption>
</figure>

<figure markdown>
  ![pgv-viewer-demo.gif](./images/pgv-viewer-demo.gif)
  <figcaption>
    Interactive HTML Viewer (<a href="./images/pgv-viewer-demo.html">Demo Page</a>)
  </figcaption>
</figure>

:sparkles: GUI (Web Application) functionality is newly added from v0.4.0

<figure markdown>
  ![pygenomeviz_gui_demo.png](https://raw.githubusercontent.com/moshi4/pyGenomeViz/main/src/pygenomeviz/gui/assets/pgv_demo.gif)
  <figcaption>
    pyGenomeViz web application example (<a href="https://pygenomeviz.streamlit.app">Demo Page</a>)
  </figcaption>
</figure>

## Installation

`Python 3.8 or later` is required for installation.

**Install PyPI package:**

    pip install pygenomeviz

**Install bioconda package:**

    conda install -c conda-forge -c bioconda pygenomeviz

**Use Docker ([Image Registry](https://github.com/moshi4/pyGenomeViz/pkgs/container/pygenomeviz)):**

Case1. Run CLI Workflow:

    docker run -it --rm ghcr.io/moshi4/pygenomeviz:latest pgv-mummer -h

Case2. Launch GUI (Web Application):

    docker run -it --rm -p 8501:8501 ghcr.io/moshi4/pygenomeviz:latest pgv-gui
