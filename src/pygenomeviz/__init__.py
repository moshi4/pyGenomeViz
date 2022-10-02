import warnings

import matplotlib as mpl

from pygenomeviz.genbank import Genbank
from pygenomeviz.genomeviz import GenomeViz
from pygenomeviz.utils import load_dataset

warnings.filterwarnings("ignore")

__version__ = "0.2.3"

__all__ = [
    "GenomeViz",
    "Genbank",
    "load_dataset",
]

# Setting matplotlib rc(runtime configuration) parameters
# https://matplotlib.org/stable/tutorials/introductory/customizing.html
mpl_rc_params = {
    # Legend
    "legend.loc": "upper left",  # Default: best
    "legend.frameon": False,  # Default: True
    "legend.handlelength": 1,  # Default: 2.0
    "legend.handleheight": 1,  # Default: 0.7
    # Savefig
    "savefig.bbox": "tight",  # Default: None
    "savefig.pad_inches": 0.5,  # Default: 0.1
    # SVG
    "svg.fonttype": "none",
}
mpl.rcParams.update(mpl_rc_params)
