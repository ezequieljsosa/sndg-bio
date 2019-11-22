FROM nanozoo/metawrap
RUN apt-get update && apt-get install -y libncurses5-dev libbz2-dev lzma liblzma-dev zlib1g-dev \
    cmake libgtk2.0-dev ncbi-blast+ libeigen3-dev unzip wget nano
RUN mkdir /app
WORKDIR /app
RUN cd /app/ && wget -O bwa-0.7.17.tar.bz2 https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2/download && tar xfv bwa-0.7.17.tar.bz2 && cd bwa-0.7.17 && make && ln -s /app/bwa-0.7.17/bwa /usr/local/bin
RUN  cd /app/ && wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && tar xfv samtools-1.9.tar.bz2 && cd samtools-1.9 && ./configure && make && make install
RUN apt-get install -y libboost-all-dev cmake
RUN wget https://bitbucket.org/berkeleylab/metabat/get/master.tar.gz && tar xzvf master.tar.gz && rm master.tar.gz &&  cd berkeleylab-metabat-* && mkdir build && cd build && cmake .. -DCMAKE_INSTALL_PREFIX=/opt/metabat && make && make test && make install
ENV PATH $PATH:/opt/metabat/bin/
WORKDIR /out

