#!/bin/bash

hmm=$1

if [ -z $1  ]
then
hmm="pfamatigrfam.hmm"
echo "Archivo no indicado. Usando pfamatigrfam.hmm"
fi

echo "Generando shorthmm para $hmm..."

grep  "NAME  "   -A15  $hmm |grep "MM    \|CONS  " -v > shorthmm
echo "--" >> shorthmm
