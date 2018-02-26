#!/bin/bash
0preprocessor.py -i $1
clean="${1/.fasta/_clean.fasta}"
#-t	Tipo de organismo (B=bacteria, A=archaea). Default: B
1rnapredictor.py -i $clean
# -r training set path
2proteinpredictor.py -i $clean
3hmmerhits.pl -i $clean
4blasthits.py -i $clean
5blasthitspostcuration.py -i $clean
6hmmerhitspostcuration.pl -i $clean
8anotator.py -i $clean -o $2 -l $3