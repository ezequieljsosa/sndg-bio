FROM ubuntu:18.04

WORKDIR /app

RUN apt-get update && apt-get install -y  openjdk-8-jdk libncurses5-dev libbz2-dev lzma liblzma-dev zlib1g-dev \
    cmake libgtk2.0-dev ncbi-blast+ libeigen3-dev unzip wget nano

RUN cd /app && wget http://cab.spbu.ru/files/release3.14.1/SPAdes-3.14.1-Linux.tar.gz && tar xfv SPAdes-3.14.1-Linux.tar.gz && rm SPAdes-3.14.1-Linux.tar.gz && ln -s /app/SPAdes-3.14.1-Linux/bin/spades.py /usr/local/bin

WORKDIR /out
CMD ["bash"]
