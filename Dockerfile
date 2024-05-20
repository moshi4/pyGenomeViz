FROM python:3.12-slim

# Install dependent aligner tools
RUN apt-get update && \
    apt-get install -y ncbi-blast+ mummer mmseqs2 progressivemauve

# Install pyGenomeViz
RUN pip install -U pip && \
    pip install pygenomeviz[gui] --no-cache-dir

# Download example dataset in advance
RUN pgv-download yersinia_phage --cache_only && \
    pgv-download mycoplasma_mycoides --cache_only

CMD ["/bin/bash"]
