FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /app
RUN apt-get update && apt-get install -y openjdk-8-jdk libncurses5-dev libbz2-dev lzma liblzma-dev zlib1g-dev \
    cmake unzip wget nano python3-pip python3 && ln /usr/bin/python3 /usr/bin/python
RUN wget -O bwa-0.7.17.tar.bz2 https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2/download \
    && tar xfv bwa-0.7.17.tar.bz2 && cd bwa-0.7.17 && make && ln /app/bwa-0.7.17/bwa /usr/local/bin && cd .. \
    && wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2 \
    && tar xfv samtools-1.12.tar.bz2 && rm samtools-1.12.tar.bz2 && cd samtools-1.12 && ./configure && make && make install && cd .. \
    && wget https://github.com/samtools/bcftools/releases/download/1.12/bcftools-1.12.tar.bz2 \
    && tar xfv bcftools-1.12.tar.bz2 && rm bcftools-1.12.tar.bz2 && cd bcftools-1.12 && ./configure && make && make install && cd .. \
    && wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz \
    && tar xfv bedtools-2.30.0.tar.gz && rm bedtools-2.30.0.tar.gz && cd bedtools2 && make && make install && cd /  \
    && rm -r /app

WORKDIR /

RUN wget -O gatk.zip "https://github.com/broadinstitute/gatk/releases/download/4.4.0.0/gatk-4.4.0.0.zip" && \
    unzip gatk.zip && rm gatk.zip && mv /gatk-4.4.0.0 /gatk
RUN pip install numpy biopython pandas tqdm
RUN apt install git
ENV MPLCONFIGDIR=/tmp

