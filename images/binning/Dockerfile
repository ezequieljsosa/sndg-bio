FROM python:3.6
WORKDIR /opt
RUN apt update && apt install -y mummer prodigal hmmer unzip curl wget python2 python2-dev python3-dev
RUN wget https://ani.jgi.doe.gov/download_files/ANIcalculator_v1.tgz && tar xfv ANIcalculator_v1.tgz && rm ANIcalculator_v1.tgz 
RUN wget https://github.com/marbl/Mash/releases/download/v1.1.1/mash-Linux64-v1.1.1.tar.gz && \
    tar xfv mash-Linux64-v1.1.1.tar.gz && cp /opt/mash-Linux64-v1.1.1/mash /usr/local/bin/ && rm -r mash-Linux64-v1.1.1.tar.gz  /opt/mash-Linux64-v1.1.1/
RUN pip3 install drep
RUN wget -O /usr/local/bin/prodigal https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux
RUN wget https://github.com/matsen/pplacer/releases/download/v1.1.alpha19/pplacer-linux-v1.1.alpha19.zip && unzip pplacer-linux-v1.1.alpha19.zip && rm pplacer-linux-v1.1.alpha19.zip 
RUN apt install -y python-pip
RUN pip2 install numpy pysam checkm-genome
RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/centrifuge/downloads/centrifuge-1.0.3-beta-Linux_x86_64.zip && unzip centrifuge-1.0.3-beta-Linux_x86_64.zip && rm centrifuge-1.0.3-beta-Linux_x86_64.zip && cd centrifuge-1.0.3-beta
RUN curl https://pyenv.run  | bash 
RUN echo 'export PATH="/root/.pyenv/bin:$PATH";eval "$(pyenv init -)";eval "$(pyenv virtualenv-init -)"' >> /root/.bashrc

RUN exec $SHELL && /root/.pyenv/bin/pyenv global 3.6.9 2.7.9

ENV PATH $PATH:/opt/ANIcalculator_v1/:/opt/pplacer-Linux-v1.1.alpha19:/opt/centrifuge-1.0.3-beta/
