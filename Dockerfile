FROM python:3.9-slim

RUN apt-get update && \
    apt-get install -y mummer mmseqs2 progressivemauve

RUN pip install -U pip && \
    pip install pygenomeviz --no-cache-dir

CMD ["/bin/bash"]
