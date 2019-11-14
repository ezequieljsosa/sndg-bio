FROM perl:5.30.0-threaded
WORKDIR /opt
RUN apt update && apt install -y nano wget git
RUN wget http://faculty.virginia.edu/wrpearson/fasta/fasta33-35/fasta-35.4.12.tar.gz && tar xfv fasta-35.4.12.tar.gz && cd fasta-35.4.12/src && make -f ../make/Makefile.linux_sse2 all
RUN git clone https://github.com/carringtonlab/TargetFinder.git
RUN apt install libssl-dev
RUN apt install  build-essential
RUN yes | cpan -i File::Which
RUN sed -i 's|#!/usr/bin/perl|#!/usr/bin/env perl|' TargetFinder/targetfinder.pl
WORKDIR /out
ENV PATH $PATH:/opt/fasta-35.4.12/bin
CMD /opt/TargetFinder/targetfinder.pl
