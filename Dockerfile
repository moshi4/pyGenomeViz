FROM python:3.10-slim

RUN apt-get update && \
    apt-get install -y mummer mmseqs2 progressivemauve

RUN pip install -U pip && \
    pip install pygenomeviz[gui] --no-cache-dir

# Download example dataset in advance
RUN pgv-download-dataset -n enterobacteria_phage

CMD ["/bin/bash"]
