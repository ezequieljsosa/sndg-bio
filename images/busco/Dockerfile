FROM ubuntu:xenial
LABEL maintainer "sujaikumar@gmail.com"


RUN apt-get update && apt-get upgrade -y -q

RUN apt-get install -y -q \
    software-properties-common \
    libboost-iostreams-dev libboost-system-dev libboost-filesystem-dev \
    zlibc gcc-multilib apt-utils zlib1g-dev python \
    cmake tcsh build-essential g++ git wget gzip perl cpanminus ncbi-blast+

COPY gm_et_linux_64.tar.gz /

RUN tar -xzf gm_et_linux_64.tar.gz

RUN wget http://bioinf.uni-greifswald.de/augustus/binaries/augustus-3.3.1.tar.gz

RUN tar -xzf augustus-3.3.1.tar.gz

RUN apt-get install bamtools -y -q

RUN wget https://github.com/Gaius-Augustus/BRAKER/archive/v2.1.2.tar.gz 2>/dev/null

RUN tar -xzf v2.1.2.tar.gz

RUN apt-get install samtools -y -q

RUN cpanm File::Spec::Functions Hash::Merge List::Util Logger::Simple \
      Module::Load::Conditional Parallel::ForkManager POSIX Scalar::Util::Numeric YAML

RUN adduser --disabled-password --gecos '' dockeruser
RUN mkdir /data
RUN chown -R dockeruser /data
RUN chmod a+w -R /augustus-3.3.1/config
USER dockeruser


COPY gm_key_64.gz /
RUN zcat /gm_key_64.gz > /home/dockeruser/.gm_key


USER root

RUN apt-get install -y python3-pip
RUN pip3 install biopython

RUN wget http://genomethreader.org/distributions/gth-1.7.1-Linux_x86_64-64bit.tar.gz
RUN tar xfv gth-1.7.1-Linux_x86_64-64bit.tar.gz



RUN wget http://mblab.wustl.edu/software/download/eval-2.2.8.tar.gz
RUN tar xfv eval-2.2.8.tar.gz


RUN rm /gm_et_linux_64.tar.gz  /v2.1.2.tar.gz  /augustus-3.3.1.tar.gz /gth-1.7.1-Linux_x86_64-64bit.tar.gz  /eval-2.2.8.tar.gz

#USER dockeruser
RUN apt install nano
RUN cpanm File::Which


RUN apt-get install -y build-essential wget git autoconf

# Install dependencies for AUGUSTUS
RUN apt-get install -y libboost-iostreams-dev zlib1g-dev
RUN apt-get install -y libgsl-dev libboost-graph-dev libsuitesparse-dev liblpsolve55-dev libsqlite3-dev libmysql++-dev
RUN apt-get install -y libbamtools-dev
RUN apt-get install -y libboost-all-dev

# Install additional dependencies for htslib and samtools
RUN apt-get install -y libbz2-dev liblzma-dev
RUN apt-get install -y libncurses5-dev

# Install additional dependencies for bam2wig
RUN apt-get install -y libssl-dev libcurl3-dev

# Build bam2wig dependencies (htslib, bfctools, samtools)
RUN git clone https://github.com/samtools/htslib.git /root/htslib
WORKDIR "/root/htslib"
RUN autoheader
RUN autoconf
RUN ./configure
RUN make
RUN make install
RUN git clone https://github.com/samtools/bcftools.git /root/bcftools
WORKDIR "/root/bcftools"
RUN autoheader
RUN autoconf
RUN ./configure
RUN make
RUN make install
RUN git clone https://github.com/samtools/samtools.git /root/samtools
WORKDIR "/root/samtools"
RUN autoheader
RUN autoconf -Wno-syntax
RUN ./configure
RUN make
RUN make install
ENV TOOLDIR="/root"

# Clone AUGUSTUS repository
WORKDIR "/root/"
RUN git clone https://github.com/Gaius-Augustus/Augustus.git augustus
#RUN wget https://github.com/Gaius-Augustus/Augustus/archive/3.3.2.tar.gz
#RUN tar xfv 3.3.2.tar.gz
#RUN mv augustus-3.3.1

# Build bam2wig
RUN mkdir -p /root/augustus/bin
WORKDIR "/root/augustus/auxprogs/bam2wig"
RUN make

# Build AUGUSTUS
WORKDIR "/root/augustus"
RUN make
RUN make install

# Test AUGUSTUS
RUN make test


USER dockeruser
WORKDIR /data

ENV BSSMDIR "/gth-1.7.1-Linux_x86_64-64bit/bin//bin/bssm"
ENV GTHDATADIR "/gth-1.7.1-Linux_x86_64-64bit/bin//bin/gthdata"
ENV PERL5LIB $PERL5LIB:/eval-2.2.8/
ENV ALIGNMENT_TOOL_PATH /gth-1.7.1-Linux_x86_64-64bit/bin/
ENV PATH $PATH:/BRAKER-2.1.2:/root/augustus/bin/:/BRAKER-2.1.2/scripts/
ENV PATH $PATH:/augustus-3.3.1/scripts:/gm_et_linux_64/gmes_petap/

ENV AUGUSTUS_BIN_PATH /root/augustus/bin/
ENV AUGUSTUS_SCRIPTS_PATH /root/augustus/scripts
ENV GENEMARK_PATH /gm_et_linux_64/
ENV BAMTOOLS_PATH /usr/bin/
ENV SAMTOOLS_PATH /usr/bin/
ENV SAMTOOLS_PATH=/usr/bin/
ENV BLAST_PATH=/usr/bin/
USER root
#RUN cp -r /root/augustus/config /home/dockeruser/augustus_config
#RUN chown dockeruser:dockeruser -R /home/dockeruser/augustus_config
RUN chown dockeruser:dockeruser -R /root

RUN apt install -y emboss
ENV PATH $PATH:/busco:/gt/current/bin


RUN wget http://genometools.org/pub/binary_distributions/gt-1.5.10-Linux_x86_64-64bit-barebone.tar.gz 
RUN tar xfv gt-1.5.10-Linux_x86_64-64bit-barebone.tar.gz 
RUN rm gt-1.5.10-Linux_x86_64-64bit-barebone.tar.gz 
RUN mv gt-1.5.10-Linux_x86_64-64bit-barebone/  /gt
WORKDIR /gt
#RUN  ln -s gt* current

RUN git clone https://gitlab.com/ezlab/busco.git  /busco
WORKDIR /busco
RUN  python setup.py install

RUN cp  /busco/config/config.ini.default  /busco/config/config.ini && chown dockeruser:dockeruser /busco/config/config.ini
RUN apt install -y hmmer
USER dockeruser
ENV AUGUSTUS_CONFIG_PATH /root/augustus/config
ADD config.ini /busco/config/config.ini
WORKDIR /out

# python /busco/scripts/run_BUSCO.py -i augustus.hints.aa --cpu 3 -o busco_search_alt  -l fungi_odb9 --mode prot
#

