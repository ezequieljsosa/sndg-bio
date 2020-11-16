FROM ubuntu:18.04

RUN apt update && apt install -y ncbi-blast+ wget nano hmmer git perl sqlite3 build-essential
RUN yes | cpan -i URI::Escape # DBI DBD::SQLite

ENV APP_NAME=Trinotate
ENV VERSION=Trinotate-v3.1.1
ENV GIT=https://github.com/Trinotate/Trinotate.git
ENV DEST=/opt/

WORKDIR $DEST

RUN git clone $GIT ; \
    cd ${APP_NAME} ; \
    git checkout tags/${VERSION} ; \
    rm -rf .git ;
RUN wget  https://github.com/trinityrnaseq/trinityrnaseq/releases/download/v2.8.6/trinityrnaseq-v2.8.6.FULL.tar.gz && tar xfv trinityrnaseq-v2.8.6.FULL.tar.gz && rm  trinityrnaseq-v2.8.6.FULL.tar.gz
RUN wget https://github.com/TransDecoder/TransDecoder/archive/TransDecoder-v5.5.0.tar.gz && tar xfv TransDecoder-v5.5.0.tar.gz && rm TransDecoder-v5.5.0.tar.gz


ENV PATH=$DEST/${APP_NAME}/:$DEST/TransDecoder-TransDecoder-v5.5.0/:$DEST/trinityrnaseq-v2.8.6/util/support_scripts/:$PATH
WORKDIR /out

#https://github.com/Trinotate/Trinotate.github.io/wiki/Software-installation-and-data-required
#https://github.com/TransDecoder/TransDecoder/wiki