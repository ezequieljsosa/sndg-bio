#!/usr/bin/python

import sys
import os
import getopt
import time
#import mysql.connector
from random import randint
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from commands import getoutput


def help():
    print  "Prediccion de rRNAs y tRNAs en una secuencia de nucleotidos.\n\
Opciones:\n\
	-i	Archivo de contigs o scaffolds. Default: contigs.fasta\n\
	-t	Tipo de organismo (B=bacteria, A=archaea). Default: B\n\
	-h 	Imprime este mensaje de ayuda\n"
	

try:
	options, args = getopt.getopt(sys.argv[1:], "i:t:h")
except getopt.GetoptError as err:
	print str(err)
	sys.exit(2)

#Parametros configurables:	
params = {}
params["i"] =  "contigs.fasta"                   	#Archivo de contigs de entrada por defecto
params["t"] =  "B" 			               	#Tipo de organismo (bacteria o archaea)



#Asigno los parametros que hayan sido definidos externamente
for option, value in options:
	if option.startswith("-"): option = option[1:]
	if option in params.keys(): params[option] = value
        if option == "h":
            help()
            sys.exit()

log = params["i"] + ".log"
if os.path.isfile(log):
	print "Ya hay un pipeline en proceso. Comenzar nuevamente (s/n)?"
	seguir = raw_input().strip()
	if seguir == 's':
		pass
	else:
		sys.exit()


now =  time.strftime("%c")

log_step = open(log, "a")
print >> log_step, "[ %s ] Comenzando Pipeline desde Paso 1. Anotacion completa.\n" % now



#Chequeo la existencia del archivo fasta
if getoutput("ls '%s' " % params["i"]) != params["i"]:
	print ("Error: No se encuentra el archivo %s \nUse el parametro -h para obtener ayuda sobre el uso del script" % params["i"])
	sys.exit(2)
if params["t"] !=  "B"  and params["t"] !=  "A":
	print ("Error: Debe ingresar un tipo correcto de organismo (A= Archaea; B= Bacteria)\nUse el parametro -h para obtener ayuda sobre el uso del script ")
	sys.exit(2)

if params["t"] == "B":
	tipornammer= "bac"
	print "Prediccion de RNAs para Bacteria"
if params["t"] == "A":
	tipornammer= "arc"
	print "Prediccion de RNAs para Archaea"
getoutput("rm " + "-r tmpbia/")
f_contigs = open(params["i"])

#Ejecucion de las predicciones
for seq_record in SeqIO.parse(f_contigs, "fasta"):
	print "	" + seq_record.id + ".fna"
	getoutput("mkdir tmpbia")
	getoutput("mkdir '%s'  " % (seq_record.id))
	fasta_file = open("%s/%s.fna" % (seq_record.id, seq_record.id), "w")
	SeqIO.write(seq_record, fasta_file, "fasta")
	fasta_file.close()
	#Deteccion de rRNAs (arc o bac)
	getoutput("rnammer -S %s  -m  lsu,ssu,tsu -g 'tmpbia/rRNA_%s.gff'   -f   '%s/rRNA_%s.fna' ' %s/%s.fna'  "  % (tipornammer, seq_record.id, seq_record.id, seq_record.id, seq_record.id, seq_record.id)  )	
	#Deteccion de tRNAs (A (archaeal) o B (bacterial))
	getoutput("tRNAscan-SE  -%s  -o 'tmpbia/tRNA_%s.report1' -f  'tmpbia/tRNA_%s.report2' '%s/%s.fna'  "   % (params["t"], seq_record.id, seq_record.id, seq_record.id, seq_record.id)  )	

	#Obtencion de archivo fasta de la salida de tRNAscan
	tRNAreport1 = open("tmpbia/tRNA_%s.report1" %  (seq_record.id) )
	tRNAs = list(tRNAreport1)
	tRNAreport1.close()

	getoutput("grep Seq 'tmpbia/tRNA_%s.report2' > 'tmpbia/tRNA_seqs_%s.report2'" %  (seq_record.id, seq_record.id) )
	tRNAreport2 = open("tmpbia/tRNA_seqs_%s.report2" %  seq_record.id  )
	tRNAseqs = list(tRNAreport2)
	tRNAreport2.close()
	
	tRNAfasta = open("%s/tRNA_%s.fna" %  (seq_record.id, seq_record.id), "w")
	i=0
	for tRNA in tRNAs[3:]:
		id, num, beg, end, type, codon, intronb, introne, score =  tRNA.strip().split("\t")[0:9]
		trnaseq_record = SeqRecord(Seq(tRNAseqs[i][5:-1]), id="tRNA_%s_%s-%s"  %  (id.strip(), beg.strip(), end.strip() ), description="/molecule=tRNA-%s /score=%s /anticodon=%s" % (type, score, codon))
		i+=1
		SeqIO.write(trnaseq_record, tRNAfasta, "fasta")
	tRNAfasta.close()
	
	# # Enmascarado de genes de RNA en el genoma
	# contigs = open(params["i"])
	# seq_iterator = SeqIO.parse(contigs, "fasta")
	# rRNAreport = open("rRNA_%s.gff" %  (file_prefix) )
	# rRNAs = list(rRNAreport)
	# rRNAreport.close()
	
	# mask_contigs = open("mask_%s.fna" % input_name, "w")
	# # mask_contigs = open("mask_%s.coords" % input_name, "w")
	
	# for record in seq_iterator:
		# seq = list(record.seq)
		# #Enmascarado de rRNA	
		# for rRNA in rRNAs[6:-1]:
			# seqname, source, feature, beg, end=  rRNA.strip().split("\t")[:5]
			# if seqname == record.id:
				# beg, end = int(beg)+50, int(end)-50
				# seq[beg-1:end] = ["N"]*(end-beg+1)
		# #Enmascarado de tRNA
		# for tRNA in tRNAs[3:]:
			# id, num, beg, end=  tRNA.strip().split("\t")[:4]
			# if id.strip() == record.id:
				# beg, end = int(beg), int(end)
				# if beg > end:
					# beg, end = end, beg
				# beg, end = beg+50, end-50
				# seq[beg-1:end] = ["N"]*(end-beg+1)
	
		# record.seq = Seq("".join(seq))
		# SeqIO.write(record, mask_contigs, "fasta")
			
	# contigs.close()
	# mask_contigs.close()
	
	#Ordenado de archivos de rRNAs y tRNAs por posicion en el genoma.
	rRNAs = open("%s/rRNA_%s.fna" % (seq_record.id, seq_record.id) )
	rRNA_list=list(SeqIO.parse(rRNAs,"fasta"))
	rRNAs.close()
	rRNA_list.sort(cmp=lambda x,y: cmp(int(x.id.strip().split("_")[-2].split("-")[0]), int(y.id.strip().split("_")[-2].split("-")[0])))
	ordered_rRNAs = open("%s/rRNA_%s.fna" % (seq_record.id, seq_record.id), "w")
	SeqIO.write(rRNA_list, ordered_rRNAs, "fasta")
	ordered_rRNAs.close()
	
	tRNAs = open("%s/tRNA_%s.fna" % (seq_record.id, seq_record.id) )
	tRNA_list=list(SeqIO.parse(tRNAs,"fasta"))
	tRNAs.close()
	tRNA_list.sort(cmp=lambda x,y: cmp(int(x.id.strip().split("_")[-1].split("-")[0]), int(y.id.strip().split("_")[-1].split("-")[0])))
	ordered_tRNAs = open("%s/tRNA_%s.fna" % (seq_record.id, seq_record.id), "w")
	SeqIO.write(tRNA_list, ordered_tRNAs, "fasta")
	ordered_tRNAs.close()

	##Almacenado en DB MySql
	#conexion = mysql.connector.connect(user='gburguener', database='BIAGenome')
	#cursor = conexion.cursor()
	#for rRNA in rRNA_list:
		# #Insertar nuevo RNA
		# add_RNA = ("INSERT INTO RNA "    "(idRNA, Tipo, Start, Stop, SecuenciaDNA, idContig) "    "VALUES (%s, %s, %s, %s, %s, %s)")
		# datos_RNA = ('Geert', 'Vanderkelen', tomorrow, 'M', date(1977, 6, 14))
		# cursor.execute(add_RNA, datos_RNA)
	#	print rRNA

	#conexion.commit()
	#cursor.close()
	#conexion.close()

f_contigs.close()

getoutput("rm "+ "-r tmpbia/")

now =  time.strftime("%c")
print >> log_step, "[ %s ] Paso 1 RNApredictor corrido sobre archivo %s para un organismo de tipo %s\n" % (now, params["i"], params["t"])
log_step.close()


