FROM ubuntu:20.04

WORKDIR /app

RUN apt-get update && DEBIAN_FRONTEND="noninteractive" apt-get install -y  openjdk-8-jdk libncurses5-dev libbz2-dev lzma liblzma-dev zlib1g-dev \
    cmake libgtk2.0-dev ncbi-blast+ libeigen3-dev unzip wget nano curl && rm -rf /var/lib/apt/lists/*

RUN  cd /app/ && wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2 \
    && tar xfv samtools-1.12.tar.bz2 && rm samtools-1.12.tar.bz2 && cd samtools-1.12 && ./configure && make && make install
RUN  cd /app/ && wget https://github.com/samtools/bcftools/releases/download/1.12/bcftools-1.12.tar.bz2 \
    && tar xfv bcftools-1.12.tar.bz2 && rm bcftools-1.12.tar.bz2 && cd bcftools-1.12 && ./configure && make && make install

RUN ln /usr/bin/python3 /usr/bin/python
RUN cd /app/ && wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz \
    && tar xfv bedtools-2.30.0.tar.gz && rm bedtools-2.30.0.tar.gz && cd bedtools2 && make && make install
RUN cd /app/ && wget -O bwa-0.7.17.tar.bz2 https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2/download \
    && tar xfv bwa-0.7.17.tar.bz2 && rm bwa-0.7.17.tar.bz2 && cd bwa-0.7.17 && make && ln -s /app/bwa-0.7.17/bwa /usr/local/bin
RUN cd /app/ && wget -O snpEff_latest_core.zip http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip/download \
    && unzip snpEff_latest_core.zip && rm snpEff_latest_core.zip
RUN cd /app/ && wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip && unzip Trimmomatic-0.38.zip && rm Trimmomatic-0.38.zip
RUN cd /app/ && wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip && \
    unzip fastqc_v0.11.7.zip && chmod 755 /app/FastQC/fastqc && ln -s /app/FastQC/fastqc /usr/local/bin && rm fastqc_v0.11.7.zip
RUN cd /app/ && wget http://cab.spbu.ru/files/release3.13.0/SPAdes-3.13.0-Linux.tar.gz && \
    tar xfv SPAdes-3.13.0-Linux.tar.gz && ln -s /app/SPAdes-3.13.0-Linux/bin/spades.py /usr/local/bin && rm SPAdes-3.13.0-Linux.tar.gz
RUN cd /app/ && wget https://github.com/ablab/quast/releases/download/quast_5.1.0rc1/quast-5.1.0rc1.tar.gz && tar xfv quast-5.1.0rc1.tar.gz \
    && rm quast-5.1.0rc1.tar.gz && cd quast-5.1.0rc1 && ./install.sh && ln -s /app/quast-5.1.0rc1/quast.py /usr/local/bin
RUN cd /app && wget https://mafft.cbrc.jp/alignment/software/mafft_7.407-1_amd64.deb && dpkg -i  mafft_7.407-1_amd64.deb && rm mafft_7.407-1_amd64.deb
RUN wget -O prinseq-lite-0.20.4.tar.gz https://sourceforge.net/projects/prinseq/files/standalone/prinseq-lite-0.20.4.tar.gz/download && \
    tar xfv prinseq-lite-0.20.4.tar.gz && rm prinseq-lite-0.20.4.tar.gz && \
    ln -s /app/prinseq-lite-0.20.4/prinseq-lite.pl /usr/local/bin/ && chmod +x /usr/local/bin/prinseq-lite.pl
RUN cd /app && wget https://github.com/broadinstitute/gatk/releases/download/4.2.0.0/gatk-4.2.0.0.zip && unzip gatk-4.2.0.0.zip \
    && rm gatk-4.2.0.0.zip && ln /app/gatk-4.2.0.0/gatk /usr/local/bin/gatk

RUN apt-get update && apt install -y python3-pip
RUN pip3 install numpy
RUN pip3 install pandas
RUN pip3 install scipy
RUN pip3 install tqdm
RUN pip3 install pyvcf 
RUN pip3 install biopython
RUN pip3 install goatools
RUN pip3 install bcbio-gff
RUN pip3 install "websockets<8"
RUN pip3 install Cython
#RUN pip3 install https://github.com/pysam-developers/pysam/archive/master.zip
RUN pip3 install hgvs

RUN apt install -y cd-hit
RUN pip3 install "matplotlib<3.1"

ENV TRIMMOMATIC /app/Trimmomatic-0.38/trimmomatic-0.38.jar
ENV PICARD /app/picard.jar
ENV SNPEFF /app/snpEff/snpEff.jar

ADD run.sh /app

RUN pip3 install sndg-bio==0.1.71
USER sndg
#ENTRYPOINT ["/app/run.sh"]
WORKDIR /out
CMD ["bash"]
