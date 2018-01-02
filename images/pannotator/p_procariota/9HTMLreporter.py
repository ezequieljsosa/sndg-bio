#!/usr/bin/python

import sys
import getopt
from random import randint
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from commands import getoutput
import re

try:
	options, args = getopt.getopt(sys.argv[1:], "i:o:l:")
except getopt.GetoptError as err:
	print str(err)
	sys.exit(2)

#Parametros configurables:	
params = {}
params["i"] = "genoma.fasta"			#Archivo fasta del genoma
params["o"] = "Bacteria"					#Organismo


#Asigno los parametros que hayan sido definidos externamente
for option, value in options:
	if option.startswith("-"): option = option[1:]
	if option in params.keys(): params[option] = value

#Chequeo la existencia de los archivos
if getoutput("ls '%s' " % params["i"]) != params["i"]:
	print ("Error: No se encuentra el archivo %s \n" % params["i"])    
	sys.exit(2)

input_name = ".".join(params["i"].split(".")[:-1])
params["o"] = "_".join(params["o"].split() )

getoutput("cp -r ../Anotator/htmlfiles %s_files/  "  % params["o"] )

#Apertura y parseado del archivo TBL
# HTML = open("%s.html" % params["o"] )
# TBL = open("ncbi_%s.tbl" % input_name)



# print >> TBL, "%s	%s	gene" % (beg, end)
# print >> TBL, "			locus_tag	%s_%s" % (params["l"], i)
# print >> TBL, "			gene	%s"  % gene
# print >> TBL, "%s	%s	CDS" % (beg, end)
# print >> TBL, "			protein_id	gnl|bia|%s_%s" % (params["l"], i)
# print >> TBL, "			product	%s" % product
# print >> TBL, "			note	%s" % dict_blast_lines_cog[ORF.id]
# print >> TBL, "			EC_number	%s" % EC

	
# #tRNAs predichos
# print "RNAs..."
# print >> TBL, "%s	%s	gene" % (beg, end)
# print >> TBL, "			locus_tag	%s_%s" % (params["l"], i)
# print >> TBL, "%s	%s	tRNA" % (beg, end)
# print >> TBL, "			product	%s" % tRNA_product


# #rRNAs predichos
# while True:
# print >> TBL, "%s	%s	gene" % (beg, end)
# print >> TBL, "			locus_tag	%s_%s" % (params["l"], i)
# print >> TBL, "%s	%s	rRNA" % (beg, end)
# print >> TBL, "			product	%s" % rRNA_product
	
# TBL.close()
