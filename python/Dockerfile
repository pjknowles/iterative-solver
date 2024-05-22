FROM continuumio/miniconda3
SHELL ["/bin/bash", "--login", "-c"]
RUN conda init bash
COPY requirements.txt .
RUN conda create -c conda-forge -n iterative-solver -y --file requirements.txt
RUN echo 'conda activate iterative-solver' >> ~/.bashrc