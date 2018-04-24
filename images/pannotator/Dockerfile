FROM ubuntu:16.04

ENV http_proxy http://proxy.fcen.uba.ar:8080
ENV https_proxy http://proxy.fcen.uba.ar:8080

RUN apt-get update && apt-get install -y \
   build-essential libpthread-stubs0-dev \
   libxml-simple-perl libperl4-corelibs-perl \
   liblog-log4perl-perl unzip wget gzip \
   python-biopython \
   ncbi-blast+ hmmer bioperl tigr-glimmer blast2 wget python-pip

WORKDIR /app

ADD ./data/ /app/

RUN cd /app/hmmer-2.3.1 && ./configure && make && \
    cd /app/ber/praze && make && \
    cd /app/tRNAscan-SE-1.3.1/ && make && make install && make clean && \
    cd /app/ELPH/sources && make

RUN pip install sndg-bio==0.1.15

ADD ./README.md /app/README.md
ADD ./p_procariota /app/p_procariota
ADD ./run.sh /app/run.sh
RUN chmod +x /app/p_procariota/annotate.sh && chmod +x /app/tbl2asn

ENV PERL5LIB=/app/ber/lib/:/app/tRNAscan-SE-1.3.1/tRNAscanSE/:/app/tRNAscan-SE-1.3.1/bin
ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/app/ber/:/app/ber/bin/:/app/ber/praze:/app/tRNAscan-SE-1.3.1/tRNAscanSE/:/app/hmmer-2.3.1/src/:/app/tRNAscan-SE-1.3.1/bin:/app/rnammer/:/app/ELPH/sources:/app/prinseq-lite-0.20.4/:/app/signalp-4.1/:/app/tmhmm-2.0c/bin:/app/LipoP1.0a/:/app/p_procariota/:/app

ENTRYPOINT ["/app/run.sh"]
CMD ["cat","/app/README.md"]

# docker build -t pannotator .
# docker tag pannotator $DOCKER_ID_USER/pannotator
# docker push $DOCKER_ID_USER/pannotator