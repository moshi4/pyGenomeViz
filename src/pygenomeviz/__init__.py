import matplotlib as mpl

from pygenomeviz.genomeviz import GenomeViz

__version__ = "1.2.0"

__all__ = [
    "GenomeViz",
]

###########################################################
# Matplotlib Runtime Configuration
###########################################################

# Setting matplotlib rc(runtime configuration) parameters
# https://matplotlib.org/stable/tutorials/introductory/customizing.html
mpl_rc_params = {
    # # Savefig
    "savefig.bbox": "tight",  # Default: None
    "savefig.pad_inches": 0.5,  # Default: 0.1
    # SVG
    "svg.fonttype": "none",
}
mpl.rcParams.update(mpl_rc_params)
