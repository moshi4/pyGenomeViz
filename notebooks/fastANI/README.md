# FastANI Visualization Example

Python script (visualize.py) that can plot the FastANI visual result figure using pyGenomeViz.

## Examples

Installation of pyGenomeViz is required to run examples.

```shell
pip install pyGenomeViz
```

It is also necessary to create FastANI visual file in advance.

```shell
fastANI -q B_quintana.fna -r B_henselae.fna --visualize -o fastani.out
```

### Example 1

colormap="hsv", link_color="grey", curve=False (Default)

```shell
python visualize.py B_quintana.fna B_henselae.fna fastani.out.visual example01.png
```

![example01.png](https://raw.githubusercontent.com/moshi4/pyGenomeViz/main/notebooks/fastANI/example01.png)  

### Example 2

colormap="viridis", link_color="red", curve=True

```shell
python visualize.py B_quintana.fna B_henselae.fna fastani.out.visual example02.png --cmap viridis --link_color red --curve 
```

![example02.png](https://raw.githubusercontent.com/moshi4/pyGenomeViz/main/notebooks/fastANI/example02.png)  
