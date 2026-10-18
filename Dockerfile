# Reproducible environment for the full pipeline, including AutoDock Vina and
# Open Babel (step 06). Build and run:
#
#   docker build -t hcc-screen .
#   docker run --rm -v "$PWD:/work" hcc-screen python scripts/run_pipeline.py
#
# Large inputs stay outside the image: mount the repository at /work and place
# the COCONUT library at data/ligands_phytochemical/ as data/README.md describes.
FROM mambaorg/micromamba:1.5.8

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

ENV PATH=/opt/conda/bin:$PATH
WORKDIR /work
CMD ["python", "scripts/run_pipeline.py", "--list"]
