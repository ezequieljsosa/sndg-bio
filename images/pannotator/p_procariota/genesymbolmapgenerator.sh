#!/bin/bash

echo "Extrayendo encabezados de trembl...";

grep "GN=" uniprot_sprot.fasta > uniprotheaders;
grep "GN=" uniprot_trembl.fasta >> uniprotheaders;

echo "Generando mapa de simbolos de genes...";
./genesymbolmapgenerator.py;

echo "Dividiendo en archivos..."
for i in {A..Z}
    do
    for j in {0..9}
    do
	echo "	Generando genesymbolmap.$i$j"
	egrep "^$i$j" genesymbolmap > genesymbolmap.$i$j
    done 

    for j in {A..Z}
    do
        echo "	Generando genesymbolmap.$i$j"
        egrep "^$i$j" genesymbolmap > genesymbolmap.$i$j
    done 

done
exit 0

