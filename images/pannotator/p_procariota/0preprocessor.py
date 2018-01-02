#!/usr/bin/python

import sys
import getopt
from random import randint
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from commands import getoutput



def help():
    print  "Preprocesamiento y chequeo de archivos a anotar.\n\
Opciones:\n\
	-i	Archivo de contigs o scaffolds. Default: contigs.fasta\n\
	-p      La secuencia es proteica. T o F Default: F\n\
	-t	Tag para reemplazar los encabezados. Default: BIA\n\
	-h 	Imprime este mensaje de ayuda\n"
	

#Parametros configurables:	
params = {}
params["i"] =  "contigs.fasta"                   	#Archivo de contigs de entrada por defecto
params["t"] =  "BIA"        		           	#Tag de encabezado por defecto
params["p"] =  "F"        		           	#Tipo de secuencia por defecto


try:
	options, args = getopt.getopt(sys.argv[1:], "i:p:t:h")
except getopt.GetoptError as err:
	print str(err)
	sys.exit(2)

#Asigno los parametros que hayan sido definidos externamente
for option, value in options:
	if option.startswith("-"): option = option[1:]
	if option in params.keys(): params[option] = value
        if option == "h":
            help()
            sys.exit()

input_name = ".".join(params["i"].split(".")[:-1])


#Chequeo la existencia del archivo fasta
if getoutput("ls '%s' " % params["i"]) != params["i"]:
	print ("Error: No se encuentra el archivo %s \nUse el parametro -h para obtener ayuda sobre el uso del script" % params["i"])    
	sys.exit(2)


aa = ""  #"-noniupac"
if params["p"] == "T":
	aa =  "-aa"


#Chequeo y reformateo de archivo fasta
print getoutput("prinseq-lite.pl %s -out_good %s_clean  -out_bad null -seq_case  upper -rm_header -seq_id  %s_  -seq_id_mappings  %s.header_mapping   -fasta  %s" % (aa, input_name, params["t"], params["i"], params["i"]))


if params["p"] == "T":
	getoutput( "mv %s_clean.fasta %s_clean.faa" % (input_name,input_name))

	


