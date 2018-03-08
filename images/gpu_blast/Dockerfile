# Use an official Python runtime as a parent image
FROM nvidia/cuda

ENV http_proxy http://proxy.fcen.uba.ar:8080
ENV https_proxy http://proxy.fcen.uba.ar:8080
ENV ftp_proxy http://proxy.fcen.uba.ar:8080

RUN apt-get update && apt-get install -y g++-4.9 gcc-4.9 wget

RUN update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-5 10 && \
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 10 && \
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-4.9 20 && \
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.9 20 && \
    update-alternatives --set gcc "/usr/bin/gcc-4.9" && \
    update-alternatives --set g++ "/usr/bin/g++-4.9"


WORKDIR /app

ADD ./data/ /app/

RUN cd /app/ && chmod +x /app/install && /app/install

RUN  ln /app/ncbi-blast-2.2.28+-src/c++/GCC493-ReleaseMT64/bin/makeblastdb /usr/local/bin/ && \
     ln /app/ncbi-blast-2.2.28+-src/c++/GCC493-ReleaseMT64/bin/blastp /usr/local/bin/


ADD ./README.md /app/README.md
ADD ./run.sh /app/run.sh

#ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/app/ber/:/app/ber/bin/:/app/ber/praze:/app/tRNAscan-SE-1.3.1/tRNAscanSE/:/app/hmmer-2.3.1/src/:/app/tRNAscan-SE-1.3.1/bin:/app/rnammer/:/app/ELPH/sources:/app/prinseq-lite-0.20.4/:/app/signalp-4.1/:/app/tmhmm-2.0c/bin:/app/LipoP1.0a/:/app/p_procariota/:/app

ENTRYPOINT ["/app/run.sh"]
CMD ["cat","/app/README.md"]
