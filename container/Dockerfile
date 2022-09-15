# base image
FROM continuumio/miniconda3
#FROM ubuntu:xenial

# install basic libraries and tools
RUN apt update && apt install -y procps wget gzip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# configure conda channels
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --add channels default

# regular conda stuff (w/ fixed tool version for minimap2)
RUN conda install -y minimap2=2.24 
RUN conda clean -a