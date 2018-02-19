#!/usr/bin/python

import re
import os
import sys
import time
import datetime
import getopt
from random import randint
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.NCBIStandalone import BlastParser
from Bio.Blast.NCBIStandalone import Iterator
from Bio.Blast import NCBIXML
from commands import getoutput
from aux_libs import tiempo_transcurrido
from aux_libs import cpu_count


def help():
    print  "Busqueda de funciones proteicas por similitud de secuencias\n\
Opciones:\n\
	-i	Archivo de contigs o scaffolds. Default: contigs.fasta\n\
	-c	Base de datos para busqueda de ortologos. Default: /data/cog/myva\n\
	-w	Base de datos para detalle de ortologos. Default: /data/cog/whog\n\
	-p	Base de datos para busqueda de perfiles enzimaticos. Default: /data/priamrpsdb/priam\n\
	-u	Base de datos para busqueda de funcion. Default: /data/uniprot/uniref50/uniref50.fasta\n\
	-s 	Base de datos secundaria para busqueda de funcion. Default: [vacio]\n\
	-e	Umbral de e-value para las busquedas. Default: 1e-5\n\
	-a	Numero de procesadores a usar. Default: maximo disponible\n\
	-h 	Imprime este mensaje de ayuda\n"

try:
	options, args = getopt.getopt(sys.argv[1:], "i:c:p:u:e:p:s:a:w:h")
except getopt.GetoptError as err:
	print str(err)
	sys.exit(2)


#Parametros configurables:	
params = {}
params["i"] = "orfs.fasta"               	#Sequencias de entrada por defecto
params["c"] = "/data/cog/myva"	#Base de datos por defecto de COG
params["w"] = "/data/cog/whog"	#Base de datos por defecto de COG
params["p"] = "/data/priamrpsdb/priam"	#Base de datos por defecto de Priam
params["u"] = "/data/uniprot/uniref50/uniref50.fasta"	#Base de datos por defecto de NCBI uniref100
params["s"] = ""									#Base de datos secundaria por defecto
params["e"] = "1e-5"                   				#Umbral de e-value
params["a"] =  "%d" % cpu_count()    					#Numero de procesadores

#Asigno los parametros que hayan sido definidos externamente
for option, value in options:
	if option.startswith("-"): option = option[1:]
	if option in params.keys(): params[option] = value
	if option == "h":
		help()
		sys.exit()

log =  params["i"] + ".log"
paso_3 = False
warning = False

if re.search("\.faa$", params["i"]):
	paso_3 = True

if not os.path.isfile(log):
	print "Advertencia. Esta iniciando el pipeline desde el Paso 5."
	warning  = True
else:
	if int(getoutput("grep \"Paso 1\" %s -c" % log)) == 0 and int(getoutput("grep \"Paso 2\" %s -c" % log)) == 0 and int(getoutput("grep \"Paso 3\" %s -c" % log)) >= 1:
		paso_3 = True	
		pass	
	else:
		if int(getoutput("grep \"Paso 1\" %s -c" % log)) == 0:
			print "Advertencia: el Paso 1 no parece haber sido ejecutado.\n"
			warning  = True
		if int(getoutput("grep \"Paso 2\" %s -c" % log)) == 0:
			print "Advertencia: el Paso 2 no parece haber sido ejecutado.\n"
			warning  = True
		if int(getoutput("grep \"Paso 3\" %s -c" % log)) == 0:
			print "Advertencia: el Paso 3 no parece haber sido ejecutado.\n"
			warning  = True
		if int(getoutput("grep \"Paso 4\" %s -c" % log)) == 0:
			print "Advertencia: el Paso 4 no parece haber sido ejecutado.\n"
			warning  = True
#	if int(getoutput("grep \"Paso 5\" %s -c" % log)) >= 1:
#		print "Advertencia: ya se ejecuto el Paso 5.\n"
if warning == True:
	print "Continuar de todos modos (s/n)?"
	respuesta = raw_input()
	if respuesta == "s":
		pass
	else:
		exit(0)

log_step = open(log, "a")


#Chequeo la existencia del archivo fasta
if not os.path.exists(params["i"]):
	print ("Error: No se encuentra el archivo %s \nUse el parametro -h para obtener ayuda sobre el uso del script" % params["i"])
	sys.exit(2)

#Parseo del archivo fasta original


print "Anotacion con bases de datos"
getoutput("mkdir tmpbia")


if paso_3 == True:
	mock_contig = params["i"].split(".")[0];
	params["i"] = mock_contig + ".fna";

f_contigs = open(params["i"])
for seq_record in SeqIO.parse(f_contigs, "fasta"):
	#Ejecucion de blastp contra COG
	queryfile = open('%s/orfs_%s.faa' % (seq_record.id, seq_record.id))
	queryrecs = list(SeqIO.parse(queryfile, "fasta"))
	queryfile.close()
	print "	Procesando " +  seq_record.id
	
	print "\tBusqueda de ortologos en COG"
	blastout = 'tmpbia/cog_%s.blast' % seq_record.id
	
	orfs_count = int(getoutput("grep '>' -c  '%s/orfs_%s.faa' " % (seq_record.id, seq_record.id) ))
	
	if os.path.exists(blastout):
		print '\tUsando archivo %s previamente calculado' % blastout
	else:
		getoutput("blastall -p blastp -i '%s/orfs_%s.faa' -d %s -e %s -a %s -v 30 -b 30 -o %s  " % (seq_record.id,seq_record.id,params["c"],params["e"],params["a"],blastout))
	

	#Parseo del archivo blast generado
	fileoutHandler =  open(blastout)
	blastparser = BlastParser()
	blast_iterator = Iterator(fileoutHandler, blastparser)
	cog_whog = open(params["w"])
	whog_lines = cog_whog.readlines()
	cog_whog.close()
	hits = []
	i=0
	for blast_record in blast_iterator:
		queryrec=queryrecs[i]
		i+=1
		for alignment in blast_record.alignments: 	#Para cada alineamiento de blastout
			sbjctID = alignment.title[1:]	#Obtengo ID de la secuencia subject
			j=0
			while j < len(whog_lines):
				if whog_lines[j][0] == "[":
					sbjctCOG = whog_lines[j].strip()
				elif re.match(".*"+sbjctID+".*", whog_lines[j] ):
					break
				j+=1
			if j==len(whog_lines):
				sbjctCOG = ""
			for hsp in alignment.hsps:                     		#Para cada hsp en cada alineamiento
				cov = int(100*(float(hsp.query_end)-float(hsp.query_start)+1.0)/len(queryrec.seq))  #Calculo la cobertura
				ident = int(float(hsp.identities[0])/float(hsp.identities[1])*100)    #Porcentaje de identidad
				score = float(hsp.score) 	#Score
				evalue = hsp.expect          	#E-value
				dfrom = int(hsp.query_start)	#Posicion de inicio y fin del hsp en el query
				dto = int(hsp.query_end)
				tfrom = int(hsp.sbjct_start)  		#Posicion de inicio y fin del hsp en el hit
				tto = int(hsp.sbjct_end)
				#Genero una lista con los datos del parseo para cada alineamiento
				hits.append((queryrec.id, queryrec.description,sbjctID,sbjctCOG, dfrom,dto,tfrom,tto,ident,cov,evalue,score))

	#Salida del ID del query, ID y cadena del pdb, posiciones de los hsp, % de identidad, cobertura y e-value
	COG = open("%s/cog_%s.blast" % (seq_record.id,seq_record.id), "w")
	print >> COG, "ID del query        Descripcion        ID del hit        COG        Hsp del query        Hsp del hit        % de identidad        Cobertura        e-value"
	for hit in hits:
		print >> COG, "\t".join([ str(f) for f in hit[1:-1]])
	
	COG.close()
	fileoutHandler.close()

	
	#Nueva ejecucion de blastp, de cada best hit de COG contra el genoma, en busca de BBHs
	print "\tBusqueda de BBH"
	fasta_cog = open("%s/cog_%s.faa" %  (seq_record.id,seq_record.id), "w")
	fasta_cog_db =  open(params["c"])
	cog_db_iterator = SeqIO.parse(fasta_cog_db, "fasta")
	blast_out_cog = open("%s/cog_%s.blast" % (seq_record.id,seq_record.id))
	blast_lines_cog = blast_out_cog.readlines()
	blast_out_cog.close()
	
	j=1
	best_blast_lines_cog =  {}
	while j < len(blast_lines_cog):
		id = blast_lines_cog[j].split("\t")[0]
		best_blast_lines_cog[ id ] = blast_lines_cog[j].split("\t")[1].strip()
		while   j < len(blast_lines_cog) and id == blast_lines_cog[j].split("\t")[0] :
			j+=1

	j=1
	for cogseq in cog_db_iterator:
		if cogseq.id in best_blast_lines_cog.values():
			SeqIO.write(cogseq, fasta_cog, "fasta")
	fasta_cog.close()
	fasta_cog_db.close()
	blastout =  'tmpbia/bbhcog_%s.blast' % seq_record.id
	
	if not os.path.exists("tempgenoma"):
		getoutput("cat  */orfs_*.faa > tempgenoma")
	
	print getoutput("formatdb -i  tempgenoma " )
	getoutput("blastall -p blastp -i '%s/cog_%s.faa' -d tempgenoma -e %s -a 1 -v 30 -b 30 -o %s" % (seq_record.id,seq_record.id,params["e"],blastout))

	
	
	
	#Parseo del archivo blast generado
	fileoutHandler =  open(blastout)
	blastparser = BlastParser()
	blast_iterator = Iterator(fileoutHandler, blastparser)
	fasta_cog_db =  open("%s/cog_%s.faa" % (seq_record.id,seq_record.id))
	cogdbrecs = list(SeqIO.parse(fasta_cog_db, "fasta"))
	BBHCOG = open("%s/bbhcog_%s.blast" % (seq_record.id, seq_record.id), "w")
	
	#Guardado del ID del query, ID y producto de COG, posiciones de los hsp, % de identidad, cobertura y e-value, para los BBH
	hits = []
	i=0
	for blast_record in blast_iterator:
		cogdbrec = cogdbrecs[i]
		i+=1
		if len(blast_record.alignments) > 0:
			bestalignment = blast_record.alignments[0] 		#El mejor alineamiento de cada record
			sbjctID = bestalignment.title[1:]	#Obtengo ID de la secuencia subject (ORF predicho)
			queryID = cogdbrec.id								#Obtengo ID de la secuencia query (miembro de COG)		
			if best_blast_lines_cog.has_key(sbjctID) and (best_blast_lines_cog[ sbjctID ] == queryID):
				BestCOGLine =  getoutput("grep '%s' '%s/cog_%s.blast'  " % (sbjctID, seq_record.id, seq_record.id)).split("\n")[0]	
				#idORF, idHit, Nombre, beg, end, = BestCOGLine.split("\t")
				print >> BBHCOG, BestCOGLine 	
	fileoutHandler.close()
	fasta_cog_db.close()
	BBHCOG.close()


	#Ejecucion de blastp contra  Uniref100 para ORFs curados
	print "\tActualizacion de evidencia de similitud con Blastp"
	orffile = open("%s/orfs_%s.faa" % (seq_record.id,seq_record.id))
	ORFs = list(SeqIO.parse(orffile, "fasta"))
	curatedorffile = open("%s/curatedorfs_%s.faa" % (seq_record.id,seq_record.id), "w")	
	for ORF in ORFs:
		if  re.search("curated", ORF.description):
			SeqIO.write(ORF, curatedorffile, "fasta")
	curatedorffile.close()
	orffile.close()

	#Parseo del archivo fasta curado
	if paso_3:
		getoutput("cp %s/orfs_%s.faa %s/curatedorfs_%s.faa" % (seq_record.id, seq_record.id, seq_record.id, seq_record.id))
	
	queryfile = open("%s/curatedorfs_%s.faa" % (seq_record.id,seq_record.id))
	curatedqueryrecs = list(SeqIO.parse(queryfile, "fasta"))
	queryfile.close()
	
	blastout = 'tmpbia/curatednr_%s.blast' % seq_record.id

	if os.path.exists(blastout):
		print '\tUsando archivo %s previamente calculado' % blastout
	else:
		getoutput("blastall -p blastp -i '%s/curatedorfs_%s.faa' -d %s -e %s -a %s -v 20 -b 20 -o %s" % (seq_record.id, seq_record.id, params["u"], params["e"], params["a"],blastout))
	#Parseo del archivo blast generado
	fileoutHandler =  open(blastout)
	blastparser = BlastParser()
	blast_iterator = Iterator(fileoutHandler, blastparser)

	hits = []
	i=0
	for blast_record in blast_iterator:
		queryrec=curatedqueryrecs[i]
		i+=1
		for alignment in blast_record.alignments: 	#Para cada alineamiento de blastout
			sbjctID = alignment.title.split(" ")[0][1:]	#Obtengo ID de la secuencia subject
			sbjttProd = " ".join(alignment.title.split("[")[0].strip().split(" ")[1:])
			for hsp in alignment.hsps:                     		#Para cada hsp en cada alineamiento
				cov = int(100*(float(hsp.query_end)-float(hsp.query_start))/len(queryrec.seq))  #Calculo la cobertura
				ident = int(float(hsp.identities[0])/float(hsp.identities[1])*100)    #Porcentaje de identidad
				score = float(hsp.score) 	#Score
				evalue = hsp.expect          	#E-value
				dfrom = int(hsp.query_start)	#Posicion de inicio y fin del hsp en el query
				dto = int(hsp.query_end)
				tfrom = int(hsp.sbjct_start)  		#Posicion de inicio y fin del hsp en el hit
				tto = int(hsp.sbjct_end)
				hlen = int(alignment.length)
				#Genero una lista con los datos del parseo para cada alineamiento
				hits.append((queryrec.id, queryrec.description,sbjctID,sbjttProd,dfrom,dto,tfrom,tto,ident,cov,evalue,hlen,score))
	
	fileoutHandler.close()

	#Salida del ID del query, ID y cadena del pdb, posiciones de los hsp, % de identidad, cobertura y e-value
	blast_out_nr = open("%s/curatednr_%s.blast" % (seq_record.id,seq_record.id), "w")
	print >> blast_out_nr, "ID del query        Descripcion        ID del hit        Producto        Hsp del query        Hsp del hit        % de identidad        Cobertura        e-value     lenHit"
	for hit in hits:
		print >> blast_out_nr,  "\t".join([ str(f) for f in hit[1:-1]])
	blast_out_nr.close()

	#Ejecucion de BER
	if not paso_3:
		file = open(blastout)
		size = file.readlines()
		file.close()
		
		if  len(size) > 0:
			getoutput("formatdb -i  '%s/orfs_%s.flanking.fna'  -o T -p F  "   %  (seq_record.id,seq_record.id) )
			"\tActualizacion de evidencia de similitud con BER"
			getoutput("ber -i %s -d '%s' -f formatdb -D  '%s/orfs_%s.flanking.fna'   -F formatdb  -m  '%s/orfs_%s.map'     -e 1e-5 -E 1e-5  -n  20 -N 0"    %   (blastout, params["u"], seq_record.id,  seq_record.id, seq_record.id, seq_record.id))
		getoutput("mv out '%s/bercurated_%s' " % (seq_record.id, seq_record.id))
			
		
		blastout = 'tmpbia/sec_%s.blast' % seq_record.id
		#Busqueda secundaria (en referencia o similar).
	if params["s"]:
		try:
			print "\tBusqueda de evidencia secundaria en %s" % params["s"]
			getoutput("formatdb -i %s -o T"  % params["s"])
			getoutput("blastall -p blastp -i '%s/orfs_%s.faa' -d %s -e %s -a %s -v 20 -b 20 -o %s" % (seq_record.id,seq_record.id,params["s"],params["e"],params["a"],blastout))
		except:
			print "Problema con archivo %s" % params["s"]
			exit(1)
		#Parseo del archivo blast generado
		fileoutHandler =  open(blastout)
		blastparser = BlastParser()
		blast_iterator = Iterator(fileoutHandler, blastparser)
		
	
		hits = []
		i=0
		for blast_record in blast_iterator:
			queryrec=queryrecs[i]
			i+=1
			for alignment in blast_record.alignments: 	#Para cada alineamiento de blastout
				sbjctID = alignment.title.split(" ")[0][1:]	#Obtengo ID de la secuencia subject
				sbjttProd = " ".join(alignment.title.split("[")[0].strip().split(" ")[1:])
				for hsp in alignment.hsps:                     		#Para cada hsp en cada alineamiento
					cov = int(100*(float(hsp.query_end)-float(hsp.query_start))/len(queryrec.seq))  #Calculo la cobertura
					ident = int(float(hsp.identities[0])/float(hsp.identities[1])*100)    #Porcentaje de identidad
					score = float(hsp.score) 	#Score
					evalue = hsp.expect          	#E-value
					dfrom = int(hsp.query_start)	#Posicion de inicio y fin del hsp en el query
					dto = int(hsp.query_end)
					tfrom = int(hsp.sbjct_start)  		#Posicion de inicio y fin del hsp en el hit
					tto = int(hsp.sbjct_end)
					hlen = int(alignment.length)
					#Genero una lista con los datos del parseo para cada alineamiento
					hits.append((queryrec.id, queryrec.description,sbjctID,sbjttProd,dfrom,dto,tfrom,tto,ident,cov,evalue,hlen,score))
		
		fileoutHandler.close()
		
		#Salida del ID del query, ID y cadena del pdb, posiciones de los hsp, % de identidad, cobertura y e-value
		blast_out_nr = open("%s/sec_%s.blast" % (seq_record.id,seq_record.id), "w")
		print >> blast_out_nr, "ID del query        Descripcion        ID del hit        Producto        Hsp del query        Hsp del hit        % de identidad        Cobertura        e-value     lenHit"
		for hit in hits:
			print >> blast_out_nr,  "\t".join([ str(f) for f in hit[1:-1]])
		blast_out_nr.close()
	

		
	#Ejecucion de rpsblast contra Priam
	print "\tBusqueda de EC number en Priam"
	#Parseo del archivo fasta original
	queryfile = open('%s/orfs_%s.faa' % (seq_record.id, seq_record.id))
	queryrecs = list(SeqIO.parse(queryfile, "fasta"))
	queryfile.close()
	
	
	rpsblastout = str("tmpbia/rps_%s.blast" % seq_record.id)
	print getoutput("rpsblast -i '%s/orfs_%s.faa' -d %s -e %s -v 20 -b 20  -m 7   -o %s" % (seq_record.id, seq_record.id,params["p"],params["e"],rpsblastout))
	#Parseo del archivo blast generado
	fileoutHandler =  open(rpsblastout)
	rpsblastparser = NCBIXML.parse(fileoutHandler)
	input_dbdir = "/".join(params["p"].split("/")[:-1])
	annotation =  input_dbdir + "/annotation_rules.xml"

	hits = []
	i=0
	
	try:
		for rpsblastrecord in rpsblastparser :
			try:
				queryrec=queryrecs[i]
			except:
				continue
			i+=1
			for alignment in rpsblastrecord.alignments: 	#Para cada alineamiento de blastout
				sbjctID = alignment.hit_id[4:]				#Obtengo ID de la secuencia subject
				ECnumber = (  getoutput("grep '%s' '%s' " % (sbjctID, annotation)  )   ).strip().split()[5][5:-1]
				for hsp in alignment.hsps:                     			#Para cada hsp en cada alineamiento
					cov = int(100*(float(hsp.query_end)-float(hsp.query_start))/len(queryrec.seq))  #Calculo la cobertura
					ident = int(float(hsp.identities)/float(hsp.align_length)*100)    				#Porcentaje de identidad
					score = float(hsp.score) 					#Score
					evalue = hsp.expect          					#E-value
					dfrom = int(hsp.query_start)			#Posicion de inicio y fin del hsp en el query
					dto = int(hsp.query_end)
					tfrom = int(hsp.sbjct_start)  				#Posicion de inicio y fin del hsp en el hit
					tto = int(hsp.sbjct_end)
					#Genero una lista con los datos del parseo para cada alineamiento
					hits.append((queryrec.id, queryrec.description,sbjctID,ECnumber, dfrom,dto,tfrom,tto,ident,cov,evalue,score))
	except:
		pass
	
	#Salida del ID del query, ID y cadena del pdb, posiciones de los hsp, % de identidad, cobertura y e-value
	EC = open("%s/ec_%s.blast" % (seq_record.id,seq_record.id), "w")
	print >> EC,  "ID del query        Descripcion        ID del hit        ECnumber        Hsp del query        Hsp del hit        % de identidad        Cobertura        e-value"
	for hit in hits:
		print >> EC,  "\t".join([ str(f) for f in hit[1:-1]])
	EC.close()	
	
	
	fileoutHandler.close()
getoutput("rm tempgenoma*")
now = time.strftime("%d/%m/%Y %H:%M:%S", time.gmtime())
print >> log_step, "[ %s ] Paso 5 BlastHitsPostCuration corrido sobre archivo %s. Base de datos para busqueda de ortologos: %s. Base de datos para busqueda de perfiles enzimaticos: %s. Base de datos para busqueda de funcion: %s. Umbral de e-value para las busquedas: %s   " % (now, params["i"], params["c"], params["p"], params["u"], params["e"] )


log_step.close()


	

