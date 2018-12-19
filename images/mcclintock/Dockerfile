FROM continuumio/miniconda

RUN apt update && apt install -y nano unzip git

WORKDIR /opt

RUN conda update conda
RUN conda config --add channels defaults
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda

RUN git clone git@github.com:bergmanlab/mcclintock.git
RUN cd mcclintock
RUN sh env_install.sh
RUN sh install.sh
RUN mkdir /out

WORKDIR /out
