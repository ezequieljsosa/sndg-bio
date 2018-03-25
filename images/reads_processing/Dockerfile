FROM ubuntu:16.04

WORKDIR /app

#ENV http_proxy http://proxy.fcen.uba.ar:8080
#ENV https_proxy http://proxy.fcen.uba.ar:8080
#ENV fps_proxy http://proxy.fcen.uba.ar:8080

RUN apt-get update && apt-get install -y  openjdk-8-jdk libncurses5-dev libbz2-dev lzma liblzma-dev zlib1g-dev \
    python-biopython cmake libgtk2.0-dev ncbi-blast+ libeigen3-dev unzip wget python-pip





ADD /data /app


RUN cd /app/ && tar xfv samtools-1.7.tar.bz2 && cd samtools-1.7 && ./configure && make && make install && \
      cd /app/ && tar xfv bcftools-1.7.tar.bz2 && cd bcftools-1.7 && ./configure && make && make install && \
      cd /app/ && tar xfv bedtools-2.27.1.tar.gz && cd bedtools2 && make && make install && \
      cd /app/ && tar xfv bwa-0.7.17.tar.bz2 && cd bwa-0.7.17 && make && ln -s /app/bwa-0.7.17/bwa /usr/local/bin && \
      cd /app/ && tar xfv GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 && \
      cd /app/ && unzip snpEff_latest_core.zip && \
      cd /app/ && tar xfv prinseq-lite-0.20.4.tar.gz  && ln -s /app/prinseq-lite-0.20.4/prinseq-lite.pl /usr/local/bin \
       && chmod +x /app/prinseq-lite-0.20.4/prinseq-lite.pl && \
      cd /app/ && unzip Trimmomatic-0.36.zip && \
      cd /app/ && unzip fastqc_v0.11.7.zip && chmod 755 /app/FastQC/fastqc && ln -s /app/FastQC/fastqc /usr/local/bin && \
      cd /app/ && tar xfv SPAdes-3.11.1-Linux.tar.gz && ln -s /app/SPAdes-3.11.1-Linux/bin/spades.py /usr/local/bin && \
      cd /app/ && tar xfv quast-4.6.3.tar.gz && cd quast-4.6.3 && ./install.sh && ln -s /app/quast-4.6.3/quast.py /usr/local/bin

RUN pip install tqdm pandas biopython goatools bcbio-gff

ENV GATK /app/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
ENV TRIMMOMATIC /app/Trimmomatic-0.36/trimmomatic-0.36.jar
ENV PICARD /app/picard.jar
ENV SNPEFF /app/snpEff/snpEff.jar

ADD run.sh /app

RUN pip install sndg-bio==0.1.14


#ENTRYPOINT ["/app/run.sh"]
CMD ["bash"]
