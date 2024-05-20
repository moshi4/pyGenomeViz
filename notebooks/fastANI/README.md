# FastANI Visualization Example

Python script (`visualize.py`) that can plot the FastANI visual result using pyGenomeViz.

## Requirements

`pyGenomeViz` installation is required to run visualize.py.

**Install PyPI package:**

    pip install pygenomeviz

**Install conda-forge package:**

    conda install -c conda-forge pygenomeviz

## Usage

    $ python visualize.py --help
    usage: python visualize.py [fasta1] [fasta2] [fastANI visual file] [plot outfile]

    Visualize conserved regions detected by fastANI

    positional arguments:
      fasta1         Input genome fasta 1
      fasta2         Input genome fasta 2
      visual_file    fastANI visual result file
      plot_outfile   Plot result outfile [*.jpg|*.png|*.svg|*.pdf]

    optional arguments:
      -h, --help     show this help message and exit
      --cmap         Colormap for matplotlib (Default: 'gist_rainbow')
      --link_color   Link color (Default: 'grey')
      --curve        Plot curved style link (Default: OFF)

## Examples

It is necessary to create FastANI visual file (*.visual) in advance.

    fastANI -q B_quintana.fna -r B_henselae.fna --visualize -o fastani.out

### Example 1

colormap=`gist_rainbow`, link_color=`grey`, curve=`False` (Default)

    python visualize.py B_quintana.fna B_henselae.fna fastani.out.visual example01.png

![example01.png](https://raw.githubusercontent.com/moshi4/pyGenomeViz/main/notebooks/fastANI/example01.png)  

### Example 2

colormap=`viridis`, link_color=`red`, curve=`True`

    python visualize.py B_quintana.fna B_henselae.fna fastani.out.visual example02.png \
                        --cmap viridis --link_color red --curve 

![example02.png](https://raw.githubusercontent.com/moshi4/pyGenomeViz/main/notebooks/fastANI/example02.png)  
