# measles-dsid-typer: zero-configuration container.
#
# Build:   docker build -t measles-dsid-typer .
# Run:     docker run --rm -v "$PWD":/data measles-dsid-typer /data/sample.fasta
#
# The base install has no third-party Python dependencies. pysam (for BAM/CRAM
# input) is installed from the python:slim image's wheels so BAM input works
# out of the box inside the container.
FROM python:3.11-slim AS base

# samtools is handy for pre-processing; pysam provides the BAM reader.
RUN apt-get update \
    && apt-get install -y --no-install-recommends samtools \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . /app

# Install with the [bam] extra so the container can read BAM/CRAM too.
RUN pip install --no-cache-dir ".[bam]"

# Default working dir for mounted data.
WORKDIR /data

ENTRYPOINT ["measles-dsid-typer"]
CMD ["--help"]
