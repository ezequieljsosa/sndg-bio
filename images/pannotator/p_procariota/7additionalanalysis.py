#!/usr/bin/python

import sys
import re
import time
import os
import getopt
from commands import getoutput
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Motif.Parsers import MAST
from random import randint
from aux_libs import cpu_count

temp_dir = "tmpbia/"
psipred_data = "/data/psipred_data/"

def help():
    print  "Predicciones de funciones, estructuras y localizaciones proteicas\n\
Opciones:\n\
	-i	Archivo de contigs o scaffolds. Default: contigs.fasta\n\
	-m	Base de datos para MAST. Default: /data/regtransbase/regtransbaselong.mast\n\
	-s	Base de datos para psipred. Default: /data/uniprot/uniref50/uniref50filt.fasta\n\
	-d	Base de datos para disopred.  Default: /data/uniprot/uniref90/uniref90.fasta\n\
	-p	Base de datos para Ps_scan. Default: /data/prosite/prosite.dat\n\
	-g	Tipo de bacteria (gram+, gram-). Default: gram+ 	\n\
	-q 	Analisis rapido (sin analisis estructural). T o F. Default: F \n\
	-h 	Imprime este mensaje de ayuda\n"
	
def disorder(query, database):
	basename = ".".join(query.split(".")[:-1])
	tempdir_basename = temp_dir + "/" + basename.split("/")[-1]
	getoutput ("rm '%s.disopred' " % basename)

	for orf in SeqIO.parse(query, "fasta"):
		tempdir_baseorf = tempdir_basename + "_" + orf.id 
		
		temp_orf = open("%s/temp.faa" % temp_dir, "w")		
		SeqIO.write(orf, temp_orf, "fasta")
		temp_orf.close()
		
		try:		
			getoutput("bia_disopred.pl  '%s/temp.faa' '%s' %d " % (temp_dir, database, cpu_count()))
			getoutput("echo %s >> '%s.disopred' " % (orf.id, basename))
			getoutput("egrep \"#\"  -v %s/temp.faa.pbdat >> '%s.disopred' " % (temp_dir, basename) )
		except:
			print "Error durante la ejecucion de Disopred"
			sys.exit(2)
			
		
		
	return basename + ".disopred"
		
		
		
def secondary_struct(query, database):	
	basename = ".".join(query.split(".")[:-1])
	tempdir_basename = temp_dir + "/" + basename.split("/")[-1]

	pn = open("%s/temp.pn" % temp_dir, "w")
	sn = open("%s/temp.sn" % temp_dir,"w")
	print >> pn, "temp.chk"
	print >> sn, "temp.faa"
	pn.close()
	sn.close()
	getoutput ("rm '%s.psipred' " % basename)
	
	for orf in SeqIO.parse(query, "fasta"):
		tempdir_baseorf = tempdir_basename + "_" + orf.id 
		temp_orf = open("%s/temp.faa" % temp_dir, "w")
		
		SeqIO.write(orf, temp_orf, "fasta")
		temp_orf.close()
		
		if not os.path.exists( '%s.chk' %  tempdir_baseorf): 		
			try:		
				getoutput("blastpgp -a %d -b 0 -j 3 -h 0.001 -d '%s' -i '%s/temp.faa' -C '%s.chk' &>'%s/temp_blastpgp' " % (cpu_count(), database, temp_dir, tempdir_baseorf, temp_dir) )
			except:
				print "Error durante la ejecucion de Psi-blast"
				sys.exit(2)
		else:
			print "			Omitiendo Psi-blast ejecutado previamente... %s " %  orf.id
	
		
		try:
			getoutput ("cp  '%s.chk'  '%s/temp.chk' " % (tempdir_baseorf, temp_dir) )			
			getoutput("makemat -P '%s/temp' " %  temp_dir)
		except:
			print "Error durante la ejecucion de makemat"
			sys.exit(2)
		
		try:
			getoutput("psipred '%s/temp.mtx' '%s/weights.dat' '%s/weights.dat2' '%s/weights.dat3' > '%s.ss' "  % (temp_dir, psipred_data, psipred_data, psipred_data, tempdir_baseorf) )
		except:
			print "Error durante la ejecucion de Psipred"

		try:
			getoutput("echo %s >> '%s.psipred' " % (orf.id, basename))
			getoutput("psipass2 '%s/weights_p2.dat' 1 1.0 1.0 '%s.ss2' '%s.ss' >> '%s.psipred' " % (psipred_data, tempdir_baseorf, tempdir_baseorf, basename) )
		except:
			print "Error durante ejecucion de Psipass2"
			sys.exit(2)
			
		
	return basename + ".psipred"


	
try:
	options, args = getopt.getopt(sys.argv[1:], "i:m:p:s:d:g:q:h")
except getopt.GetoptError as err:
	print str(err)
	sys.exit(2)
	
paso_3 = False
warning = False


#Parametros configurables:	
params = {}
params["i"] = "contigs.fasta"        							#Sequencia de entrada por defecto
params["m"] = "/data/regtransbase/regtransbaselong.mast" 			#Base de datos para MAST
params["p"] = "/data/prosite/prosite.dat" 						#Base de datos para Ps_scan
params["s"] = "/data/uniprot/uniref50/uniref50filt.fasta"			#Base de datos para psipred
params["d"] = "/data/uniprot/uniref90/uniref90.fasta"				#Base de datos para disopred
params["q"] = "F"		               						#Analisis rapido (sin estructural)
params["g"] = "gram+" 		#Tipo de bacteria

#Chequeo de parametro
if params["g"] != "gram+" and params["g"] != "gram-":
	print ("Error: Parametro incorrecto. Ingrese \"gram+\" o \"gram-\" para indicar el tipo de bacteria\nUse el parametro -h para obtener ayuda sobre el uso del script")
	sys.exit(2)

#Asigno los parametros que hayan sido definidos externamente
for option, value in options:
	if option.startswith("-"): option = option[1:]
	if option in params.keys(): params[option] = value
	if option == "h":
		help()
		sys.exit()
		
log =  params["i"] + ".log"
log_step = open(log, "a")
if re.search("\.faa$", params["i"]):
	paso_3 = True

if not os.path.isfile(log):
	print "Advertencia. Esta iniciando el pipeline desde el Paso 7."
	warning = True
else:
	if int(getoutput("grep \"Paso 1\" %s -c" % log)) == 0 and int(getoutput("grep \"Paso 2\" %s -c" % log)) == 0 and int(getoutput("grep \"Paso 3\" %s -c" % log)) >= 1:
		paso_3 = True		
		pass	
	else:
		if int(getoutput("grep \"Paso 1\" %s -c" % log)) == 0:
			print "Advertencia: el Paso 1 no parece haber sido ejecutado.\n"
			warning = True
		if int(getoutput("grep \"Paso 2\" %s -c" % log)) == 0:
			print "Advertencia: el Paso 2 no parece haber sido ejecutado.\n"
			warning = True
		if int(getoutput("grep \"Paso 3\" %s -c" % log)) == 0:
			print "Advertencia: el Paso 3 no parece haber sido ejecutado.\n"
			warning = True
		if int(getoutput("grep \"Paso 4\" %s -c" % log)) == 0:
			print "Advertencia: el Paso 4 no parece haber sido ejecutado.\n"
			warning = True
	if int(getoutput("grep \"Paso 5\" %s -c" % log)) == 0:
			print "Advertencia: el Paso 5 no parece haber sido ejecutado.\n"
			warning = True
	if int(getoutput("grep \"Paso 6\" %s -c" % log)) == 0:
			print "Advertencia: el Paso 6 no parece haber sido ejecutado.\n"
			warning = True
#	if int(getoutput("grep \"Paso 7\" %s -c" % log)) >= 1:
#		print "Advertencia: ya se ejecuto el Paso 7.\n"
if warning == True:
	print "Continuar de todos modos (s/n)?"
	respuesta = raw_input()
	if respuesta == "s":
		pass
	else:
		exit(0)

log_step = open(log, "a")


#Chequeo la existencia de los archivos
if getoutput("ls '%s' " % params["i"]) != params["i"]:
	print ("Error: No se encuentra el archivo %s \nUse el parametro -h para obtener ayuda sobre el uso del script" % params["i"])
	sys.exit(2)		

# #Ejecucion de MAST. Busqueda de motivos conocidos
# print "-Busqueda de motivos conocidos"
# getoutput("mast  %s %s " % (params["m"], params["u"]))


print "Prediccion de motivos, seniales, localizacion y estructura secundaria"

if paso_3 == True:
	mock_contig = params["i"].split(".")[0];
	params["i"] = mock_contig + ".fna";

f_contigs = open(params["i"])
for seq_record in SeqIO.parse(f_contigs, "fasta"):
	seq_file = "%s/orfs_%s.faa"  %  (seq_record.id, seq_record.id)
	
	#SignalP. Prediccion de peptidos senial reconocidos por Spase I (exportacion de proteinas).
	print "Procesando",  seq_record.id
	print "	Prediccion de peptidos senial"
	getoutput("signalp -t %s -f short %s > %s/orfs_%s.signalplong"  % (params["g"], seq_file, seq_record.id, seq_record.id) )
	getoutput("egrep  \"predictions|Y\" %s/orfs_%s.signalplong > %s/orfs_%s.signalp" % (seq_record.id, seq_record.id, seq_record.id, seq_record.id))

	#LipoP. Prediccion de peptidos senial reconocidos por Spase I  (proteinas) y  Spase II (lipoproteinas) (o localizacion citoplasmatica)
	print "	Prediccion de localizacion de proteinas"
	getoutput("LipoP  -short  %s  > %s/orfs_%s.lipoplong" % (seq_file, seq_record.id, seq_record.id) )
	getoutput("egrep  -v \"CYT|TMH\" %s/orfs_%s.lipoplong > %s/orfs_%s.lipop" % (seq_record.id, seq_record.id, seq_record.id, seq_record.id))

	#TMhmm. Prediccion de dominios extracelulares, intracelulares y transmembrana.
	print "	Prediccion de dominios extracelulares, intracelulares y transmembrana"
	getoutput("tmhmm --short  %s   >%s/orfs_%s.tmhmmlong" % (seq_file, seq_record.id, seq_record.id))
	getoutput("egrep  \"PredHel=[1-9]\" %s/orfs_%s.tmhmmlong > %s/orfs_%s.tmhmm" % (seq_record.id, seq_record.id, seq_record.id, seq_record.id))
	getoutput("rm -r TMHMM_*  &&  rm %s/*long" %  seq_record.id)
	
	if	params["q"] != "T":
		#Prediccion de estructura secundaria
		print "	Prediccion de estructura secundaria"
		psipred_file = secondary_struct(seq_file, params["s"])
		
		#Prediccion de desorden estructural
		print "	Prediccion de desorden intrinseco"
		disopred_file = disorder(seq_file, params["d"])
	
	# #Ps-scan. Busqueda de patrones de aminoacidos de diferentes funciones.
	# print "-Busqueda de patrones de aminoacidos"
	# getoutput("ps_scan -d %s %s > %s.psscan"   % (params["p"], params["i"], input_name) )

	# Prediccion de estructuras secundarias
	

now = time.strftime("%d/%m/%Y %H:%M:%S", time.gmtime())
print >> log_step, "[ %s ] Paso 7 AdditionalAnalysis corrido sobre archivo %s para bacterias %s " % (now, params["i"], params["g"])

log_step.close()
