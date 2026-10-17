# Pinned by digest so that a CI result only changes when this repository does.
# This is the image that docker.io/continuumio/miniconda3:latest resolved to on
# 2026-10-16 (built 2026-08-28), so pinning it changes nothing about what CI runs
# today. Bump it deliberately; note that the newest version tag (v25.11.1-1) is a
# different image from latest.
FROM continuumio/miniconda3@sha256:eca594d684f495c1a02beff33a9fab53aec8c5830eaf431bb149912dc6c9e4c1
SHELL ["/bin/bash", "--login", "-c"]
RUN conda init bash
COPY requirements.txt .
RUN conda create -c conda-forge -n iterative-solver -y --file requirements.txt
RUN echo 'conda activate iterative-solver' >> ~/.bashrc