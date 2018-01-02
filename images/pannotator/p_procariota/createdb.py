#!/usr/bin/python

/data/uniprot/uniref90/uniref90.fasta
ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz

/data/uniprot/gsm/uniprot_sprot.fasta
ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
/data/uniprot/gsm/uniprot_trembl.fasta
ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
/data/uniprot/gsm/genesymbolmapgenerator.py
scripts/genesymbolmapgenerator.py
/data/uniprot/gsm/genesymbolmapgenerator.sh 
scripts/genesymbolmapgenerator.sh  -> gsm/*

/data/uniprot/goa/goa_uniprot_all.gpa
ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gpa.gz
goagenerator.py  -> goa/*
scripts/goagenerator.py  

/data/cog/whog
ftp://ftp.ncbi.nih.gov/pub/COG/COG/whog

/data/cog/myva
ftp://ftp.ncbi.nih.gov/pub/COG/COG/myva
formatdb -i /data/cog/myva -o T

ls /data/databases/ec/PRIAM_MAR15/PROFILES/*.chk > priam
/data/databases/ec/PRIAM_MAR15
formatrpsdb -i /data/databases/ec/PRIAM_MAR15/priam -o T

ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
/data/pfamtigrfam/pfamatigrfam.hmm
hmmconvert TIGRFAMs_15.0_HMM.LIB > tirgfam.hmm
cat tirgfam.hmm Pfam-A.hmm > pfamatigrfam.hmm
hmmpress /data/pfamtigrfam/pfamatigrfam.hmm

/data/pfamtigrfam/INFO
wget ftp://ftp.jcvi.org/pub/data/TIGRFAMs/TIGRFAMs_15.0_INFO.tar.gz

/data/pfamtigrfam/pfam2go
http://geneontology.org/external2go/pfam2go
/data/pfamtigrfam/TIGRFAMS_GO_LINK
http://geneontology.org/external2go/tigrfams2go
/data/pfamtigrfam/shorthmm_generator.sh
scripts/shorthmm_generator.sh -> /data/pfamtigrfam/shorthmm


#Armado base

mkdir cog
mkdir pfamtigrfam
mkdir priamrpsdb
mkdir uniprot
mkdir uniprot/goa
mkdir uniprot/gsm
mkdir uniprot/uniref

cd pfamtigrfam
mkdir INFO
cd INFO
wget wget ftp://ftp.jcvi.org/pub/data/TIGRFAMs/TIGRFAMs_15.0_INFO.tar.gz
tar xfv TIGRFAMs_15.0_INFO.tar.gz
cd ..
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
gunzip Pfam-A.hmm.gz
wget ftp://ftp.jcvi.org/pub/data/TIGRFAMs/TIGRFAMs_15.0_HMM.LIB.gz
gunzip TIGRFAMs_15.0_HMM.LIB.gz
hmmconvert TIGRFAMs_15.0_HMM.LIB > tirgfam.hmm
cat tirgfam.hmm Pfam-A.hmm > pfamatigrfam.hmm
shorthmm_generator.sh
wget ftp://ftp.jcvi.org/pub/data/TIGRFAMs/TIGRFAMS_GO_LINK
wget http://www.geneontology.org/external2go/pfam2go
mv pfam2go pfam2go.txt
cd ..


cd uniprot/uniref
ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
gunzip uniref90.fasta.gz
cd ../goa
wget ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gpa.gz
gunzip goa_uniprot_all.gpa.gz
./goagenerator.py  
cd ../gsm
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
gunzip uniprot_trembl.fasta.gz
./genesymbolmapgenerator.sh
cd ..



for x in `find /home/eze/uniprot/goa/ -printf "%f\n" | grep 'gp_association.goa_uniprot.'`; do ln /home/eze/uniprot/goa/$x /home/eze/annotator/uniprot/goa/$x  ; done


# POST RSYNC
cd /data/priamrpsdb
formatrpsdb -i /data/priamrpsdb/priam -o T
hmmpress /data/pfamtigrfam/pfamatigrfam.hmm
formatdb -i /data/cog/myva -o T
makeblastdb -type prot -in /data/uniprot/uniref/uniref90.fasta







