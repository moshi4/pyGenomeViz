# pyGenomeViz

![Python3](https://img.shields.io/badge/Language-Python3-steelblue)
![OS](https://img.shields.io/badge/OS-_Windows_|_Mac_|_Linux-steelblue)
![License](https://img.shields.io/badge/License-MIT-steelblue)
[![Latest PyPI version](https://img.shields.io/pypi/v/pygenomeviz.svg)](https://pypi.python.org/pypi/pygenomeviz)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/pygenomeviz.svg?color=green)](https://anaconda.org/bioconda/pygenomeviz)  
[![CI](https://github.com/moshi4/pygenomeviz/workflows/CI/badge.svg)](https://github.com/moshi4/pygenomeviz/actions/workflows/ci.yml)

## Overview

pyGenomeViz is a genome visualization python package for comparative genomics.
It is implemented based on matplotlib, the most popular visualization library in python,
and can easily and beautifully plot genomic features and comparison results.

## Installation

**Install PyPI package:**

    pip install pygenomeviz

**Install bioconda package:**

    conda install -c conda-forge -c bioconda pygenomeviz

## Usage

### Basic Usage

The example codes shown here are also available from [jupyter notebook](https://github.com/moshi4/pyGenomeViz/blob/main/example/tutorial.ipynb).

#### Single Genome Track Visualization

```python
from pygenomeviz import GenomeViz

name, genome_size = "Tutorial 01", 5000
cds_list = ((100, 900, -1), (1100, 1300, 1), (1350, 1500, 1), (1520, 1700, 1), (1900, 2200, -1), (2500, 2700, 1), (2700, 2800, -1), (2850, 3000, -1), (3100, 3500, 1), (3600, 3800, -1), (3900, 4200, -1), (4300, 4700, -1), (4800, 4850, 1))

gv = GenomeViz()
track = gv.add_feature_track(name, genome_size)
for idx, cds in enumerate(cds_list, 1):
    start, end, strand = cds
    track.add_feature(start, end, strand, label=f"CDS{idx:02d}")

fig = gv.plotfig(dpi=100)
```

#### Multiple Genome Track & Link Visualization

```python
from pygenomeviz import GenomeViz

genome_list = (
    {"name": "genome 01", "size": 1000, "cds_list": ((150, 300, 1), (500, 700, -1), (750, 950, 1))},
    {"name": "genome 02", "size": 1300, "cds_list": ((50, 200, 1), (350, 450, 1), (700, 900, -1), (950, 1150, -1))},
    {"name": "genome 03", "size": 1200, "cds_list": ((150, 300, 1), (350, 450, -1), (500, 700, -1), (701, 900, -1))},
)

gv = GenomeViz(tick_style="axis")
for genome in genome_list:
    name, size, cds_list = genome["name"], genome["size"], genome["cds_list"]
    track = gv.add_feature_track(name, size)
    for idx, cds in enumerate(cds_list, 1):
        start, end, strand = cds
        track.add_feature(start, end, strand, label=f"gene{idx:02d}", linewidth=1, labelrotation=0, labelvpos="top", labelhpos="center", labelha="center")

# Add links between "genome 01" and "genome 02"
gv.add_link(("genome 01", 150, 300), ("genome 02", 50, 200))
gv.add_link(("genome 01", 700, 500), ("genome 02", 900, 700))
gv.add_link(("genome 01", 750, 950), ("genome 02", 1150, 950))
# Add links between "genome 02" and "genome 03"
gv.add_link(("genome 02", 50, 200), ("genome 03", 150, 300), normal_color="skyblue", inverted_color="lime")
gv.add_link(("genome 02", 350, 450), ("genome 03", 450, 350), normal_color="skyblue", inverted_color="lime")
gv.add_link(("genome 02", 900, 700), ("genome 03", 700, 500), normal_color="skyblue", inverted_color="lime")
gv.add_link(("genome 03", 900, 701), ("genome 02", 1150, 950), normal_color="skyblue", inverted_color="lime")

fig = gv.plotfig(dpi=100)
```

### Practical Usage

The example codes shown here are also available from [jupyter notebook](https://github.com/moshi4/pyGenomeViz/blob/main/example/tutorial.ipynb).

#### Single Genome Track Visualization from Genbank file

```python
from pygenomeviz import Genbank, GenomeViz, load_dataset
```

#### Multiple Genome Track & Link Visualization from Genbank files

```python
from pygenomeviz import Genbank, GenomeViz, load_dataset
```

### Customization Tips

Since pyGenomeViz is implemented based on matplotlib, users can easily customize
the figure in the manner of matplotlib. Here are some tips for figure customization.

- Add `GC content` & `GC skew` subtrack
- Add annotation (Fill Box, ROI)
- Add colorbar (Experimetal implementation)
