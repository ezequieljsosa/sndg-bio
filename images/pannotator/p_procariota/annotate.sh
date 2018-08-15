#!/bin/bash

cd /out

if [ "$1" = "" ]; then
       echo "must specify fasta file"
       exit 1
fi

if [ ! -e $1 ]; then
       echo "$1 does no exists"
       exit 1
fi

if [ "$2" = "" ]; then
       echo "must specify organism name"
       exit 1
fi
if [ "$3" = "" ]; then
       echo "must specify locus tag"
       exit 1
fi

org_type=${4:-B } # (B=bacteria, A=archaea). Default: B
gram=${5:-gram+ } # gram+ gram-

reffasta=""
if [ "$6" != "" ]; then
      reffasta="-r $6"
      if [ ! -e $6 ]; then
             echo "$6 does no exists"
             exit 1
      fi
fi

0preprocessor.py -i $1

clean="${1/.fasta/_clean.fasta}"

1rnapredictor.py -i $clean -t $org_type

2proteinpredictor.py -i $clean  $reffasta
3hmmerhits.pl -i $clean
4blasthits.py -i $clean
5blasthitspostcuration.py -i $clean
6hmmerhitspostcuration.pl -i $clean
7additionalanalysis.py -q T -g $gram -i $clean
8anotator.py -i $clean -o $2 -l $3