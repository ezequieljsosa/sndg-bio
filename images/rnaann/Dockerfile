FROM ubuntu:18.04
WORKDIR /tmp
RUN apt-get update && \
    apt-get -y install wget gcc make && \
    apt-get clean

WORKDIR /app

RUN wget http://trna.ucsc.edu/software/trnascan-se-2.0.3.tar.gz
RUN tar xfv trnascan-se-2.0.3.tar.gz && rm trnascan-se-2.0.3.tar.gz
RUN cd tRNAscan-SE-2.0 && ./configure && make && make install && make clean
RUN apt-get install -y cpanminus
RUN cpanm -i File::Copy

RUN wget eddylab.org/infernal/infernal-1.1.2.tar.gz && tar xf infernal-1.1.2.tar.gz && rm infernal-1.1.2.tar.gz
RUN cd infernal-1.1.2 && ./configure && make && make install && make clean
RUN cd infernal-1.1.2/easel && make && make install && make clean

RUN wget http://lowelab.ucsc.edu/software/snoscan.tar.gz && tar xfv snoscan.tar.gz && rm snoscan.tar.gz
RUN wget http://eddylab.org/software/squid/squid.tar.gz && tar xfv squid.tar.gz && rm squid.tar.gz
RUN cd squid-1.9g/ && ./configure && make && make install && make clean
RUN mkdir -p /root/bin/ && cd snoscan-0.9.1/  && make && make install && make clean && mv /root/bin/* /usr/local/bin/
RUN apt install -y ncbi-blast+
RUN mkdir /data && echo "#!/bin/sh\ncd /data\n" \
         "wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz;gunzip Rfam.cm.gz\n" \
         "wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin\n" \
         "wget ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.dat.gz;gunzip miRNA.dat.gz\n" \
         "wget ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz;zcat hairpin.fa.gz > miRNA.fasta\n" \
         "wget ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz;zcat mature.fa.gz >> miRNA.fasta\n" \
         "rm mature.dat.gz ;rm hairpin.dat.gz\n" \
         "cmpress Rfam.cm;makeblastdb -dbtype nucl -in miRNA.fasta" > /usr/local/bin/downloadDBs.sh && chmod +x /usr/local/bin/downloadDBs.sh

WORKDIR /out