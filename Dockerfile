FROM python:3.10-slim

RUN apt-get update && \
    apt-get install -y mummer mmseqs2 progressivemauve

RUN pip install -U pip && \
    pip install pygenomeviz[all] --no-cache-dir

CMD ["/bin/bash"]
