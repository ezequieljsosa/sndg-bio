#!/bin/bash


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

org_type=${4:-\-t B } # (B=bacteria, A=archaea). Default: B

reffasta = ""
if [ "$5" = "" ]; then
      reffasta = "-r $5"
fi

0preprocessor.py -i $1

clean="${1/.fasta/_clean.fasta}"

1rnapredictor.py -i $clean

2proteinpredictor.py -i $clean $org_type $reffasta
3hmmerhits.pl -i $clean
4blasthits.py -i $clean
5blasthitspostcuration.py -i $clean
6hmmerhitspostcuration.pl -i $clean
8anotator.py -i $clean -o $2 -l $3