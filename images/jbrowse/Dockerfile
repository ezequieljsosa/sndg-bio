# JBrowse
# VERSION 1.0
FROM node:lts
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get -qq update --fix-missing
RUN apt-get --no-install-recommends -y install git build-essential zlib1g-dev libxml2-dev libexpat-dev postgresql-client libpq-dev

RUN mkdir -p /jbrowse/ && GIT_SSL_NO_VERIFY=1 git clone --recursive https://github.com/GMOD/jbrowse.git /jbrowse/ && \
    cd /jbrowse/ && \
    git checkout 1.16.8-release

WORKDIR /jbrowse/
RUN JBROWSE_PUBLIC_PATH=/static/jbrowse/dist/ ./setup.sh 
RUN ./bin/cpanm --force JSON Hash::Merge PerlIO::gzip Devel::Size \
    Heap::Simple Heap::Simple::XS List::MoreUtils Exception::Class Test::Warn Bio::Perl \
    Bio::DB::SeqFeature::Store File::Next Bio::DB::Das::Chado && \
    rm -rf /root/.cpan/

RUN perl Makefile.PL && make && make install

