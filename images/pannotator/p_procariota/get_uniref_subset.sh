#!/bin/sh

grep $1 -i  $2 > $2.$1
TblToFasta2  $2.$1 > $2.fasta.$1
rm $2.$1
