#!/usr/bin/python

import sys
import re
import getopt
from random import randint
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Blast.NCBIStandalone import BlastParser
from Bio.Blast.NCBIStandalone import Iterator
from commands import getoutput

def help():
    print  "Prediccion de genes en una secuencia de nucleotidos.\n\
Opciones:\n\
	-i	Archivo de contigs o scaffolds. Default: contigs.fasta\n\
	-t	Umbral de entropia para ORFs detectados por long-orfs. Default: 1.15\n\
	-o	Maximo overlap entre genes. Default: 50\n\
	-l	Minima longitud de gen. Default: 110\n\
	-s	Minimo score para ser considerado gen. Default: 30\n\
	-h 	Imprime este mensaje de ayuda\n"

try:
    options, args = getopt.getopt(sys.argv[1:], "i:t:o:l:s:h")
except getopt.GetoptError as err:
    print str(err)
    sys.exit(2)

#Parametros configurables:	
params = {}
params["i"] =  "contigs.fasta"                   	#Archivo de contigs de entrada por defecto
params["t"] =  "1.15"                                      	#Umbral de entropia para ORFs detectados por long-orfs
params["o"] =  "50"						#Maximo overlap entre genes
params["l"] =  "110"						#Minima longitud de gen
params["s"] =  "30"						#Minimo score para ser considerado gen

#Asigno los parametros que hayan sido definidos externamente
for option, value in options:
    if option.startswith("-"): option = option[1:]
    if option in params.keys(): params[option] = value
    if option == "h":
        help()
        sys.exit()
    
#Chequeo la existencia del archivo fasta
if getoutput("ls '%s' " % params["i"]) != params["i"]:
    print ("Error: No se encuentra el archivo %s \nUse el parametro -h para obtener ayuda sobre el uso del script" % params["i"])    
    sys.exit(2)

file_prefix = str(randint(100000,1000000))

#Primera ronda de prediccion
print "-Prediccion de genes"
getoutput("mkdir tmpbia")
f_contigs = open(params["i"])
for seq_record in SeqIO.parse(f_contigs, "fasta"):
	getoutput("mkdir '%s' " % (seq_record.id))
	fasta_file = open("%s/%s.fna" % (seq_record.id, seq_record.id), "w")
	SeqIO.write(seq_record, fasta_file, "fasta")
	fasta_file.close()
	#Deteccion inicial de ORFs largos no solapantes para entrenamiento con el propio genoma
	getoutput("tigr-glimmer long-orfs  -n  -t %s     '%s/%s.fna'   tmpbia/%s.longorfs" % (params["t"], seq_record.id, seq_record.id, file_prefix))	
	#Extraccion y guardado de ORFs de entrenamiento de cada archivo fasta
	getoutput("tigr-glimmer extract   -t   '%s/%s.fna'  tmpbia/%s.longorfs >> tmpbia/%s.train" % (seq_record.id, seq_record.id, file_prefix, file_prefix))		
f_contigs.close()

#Entrenamiento. Construccion del Interpolated Context Model
getoutput("tigr-glimmer build-icm -r tmpbia/%s.icm < tmpbia/%s.train" % (file_prefix, file_prefix))		

totalorfs = 0
totalstarts = [0, 0, 0]
f_contigs = open(params["i"])
for seq_record in SeqIO.parse(f_contigs, "fasta"):
	#Prediccion preliminar de genes con glimmer3
	getoutput("tigr-glimmer glimmer3  -z 11 -o %s -g  %s -t %s -l '%s/%s.fna' tmpbia/%s.icm tmpbia/%s" % (params["o"], params["l"], params["s"], seq_record.id, seq_record.id, file_prefix, file_prefix))		
	getoutput("tail -n +2  tmpbia/%s.predict >  tmpbia/%s.coords"   %   (file_prefix, file_prefix))										#Obtencion de las coordenadas de la primera prediccion 
	#Obtencion de regiones upstream de los codones start
	getoutput("upstream-coords.awk  25  0  tmpbia/%s.coords | tigr-glimmer extract   '%s/%s.fna' - >> tmpbia/%s.upstream.fna"  % ( file_prefix, seq_record.id, seq_record.id, file_prefix))	
	#Calculo de proporcion de uso de codones	start
	atg, gtg, ttg = getoutput("tigr-glimmer start-codon-distrib -3 '%s/%s.fna' tmpbia/%s.coords"  %  (seq_record.id, seq_record.id, file_prefix)).split(",")				
	orfcount = float(getoutput("grep -c orf  tmpbia/%s.coords"  %  (file_prefix)) )
	totalstarts = totalstarts[0]+round(float(atg)*orfcount), totalstarts[1]+round(float(gtg)*orfcount), totalstarts[2]+round(float(ttg)*orfcount) 
	totalorfs+=orfcount	
f_contigs.close()
	
totalstarts =  (totalstarts[0])/float(totalorfs)*1000.0, (totalstarts[1])/float(totalorfs)*1000.0, (totalstarts[2])/float(totalorfs)*1000.0 
startuse = "%f,%f,%f"  % (float( round(totalstarts[0]))/1000.0, float(round(totalstarts[1]))/1000.0, float(round(totalstarts[2]))/1000.0)
flank = 300
flank_startuse =  open("flank_startuse", "w")
print >> flank_startuse, flank
print >> flank_startuse, startuse
flank_startuse.close()
 
#Calculo de perfiles de RBS (PWM)  de las regiones upstream
getoutput("elph  tmpbia/%s.upstream.fna LEN=6 |   get-motif-counts.awk >  tmpbia/%s.motif" % ( file_prefix, file_prefix) )	

#Segunda ronda de prediccion
f_contigs = open(params["i"])
for seq_record in SeqIO.parse(f_contigs, "fasta"):
	print "	Procesando",  seq_record.id
	#Prediccion refinada con el perfil RBS y el uso de codones start
	getoutput("tigr-glimmer glimmer3  -z 11 -o %s -g %s -t %s -l -b tmpbia/%s.motif  -P %s '%s/%s.fna' tmpbia/%s.icm tmpbia/%sfinal  "  % (params["o"], params["l"], params["s"], file_prefix, startuse, seq_record.id, seq_record.id, file_prefix, file_prefix))	
	#Extraccion final de ORFs
	getoutput("tail -n +2  tmpbia/%sfinal.predict > tmpbia/%sfinal.coords"    %   (file_prefix, file_prefix))	
	getoutput("tigr-glimmer extract  -w  '%s/%s.fna'  tmpbia/%sfinal.coords >  tmpbia/%sfinal.fna"  %  (seq_record.id, seq_record.id, file_prefix, file_prefix ))
	fasta_file = open("tmpbia/%sfinal.fna" % file_prefix)
	orfs_seqs = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
	fasta_file.close()
	predictfinal_file = open("tmpbia/%sfinal.predict" % file_prefix)
	orffile = open("%s/orfs_%s.faa" % (seq_record.id, seq_record.id), "w")
	mapprotseq =  open("%s/orfs_%s.map" % (seq_record.id, seq_record.id),"w")		#mapa de correspondencia entre identificadores de proteinas y nucleotidos, para uso con BER
	for line in predictfinal_file:
		if line[0] == ">":
			contigID = line[1:].strip().split()[0]
			continue
		orfID, beg, end, frame, score= filter(lambda s: len(s)>0, line.strip().split(" "))
		dnaseq = Seq(str(orfs_seqs[orfID].seq), generic_dna) 
		protseq = str(dnaseq.translate(table=11))  #'M' + str(dnaseq.translate(table=11))[1:]
		if protseq.count("X") <= 25:
			orfID = orfID.split("_")[0]
			predict=SeqRecord(Seq(protseq),id=orfID,description="%s %d %d %.2f"%(contigID, int(beg), int(end), float(score)))	
			SeqIO.write(predict, orffile, "fasta")
			print >> mapprotseq, orfID + "\t" + orfID
	orffile.close()
	mapprotseq.close()
	predictfinal_file.close()
	#Extraccion de secuencia de DNA de cada gen, flanqueada por un excedente de 300bp upstream y downstream,
	#para posterior busqueda de frameshifts y mutaciones puntuales con BER
	flankingcoords_file  = open("tmpbia/%s.flankingcoords" %  (seq_record.id), "w")
	predictfinal_file = open("tmpbia/%sfinal.predict" % file_prefix)
	for line in predictfinal_file:
		if line[0] == ">":
			continue
		orfID, beg, end, frame, score = filter(lambda s: len(s)>0, line.strip().split(" "))
		direction = int(frame)/abs(int(frame))
		beg = int(beg) - direction*flank
		end = int(end) + direction*flank
		# if beg<1:
			# beg = 1
		# elif beg > len(seq_record.seq):
			# beg = len(seq_record.seq)
		# if end<1:
			# end = 1
		# elif end > len(seq_record.seq):
			# end = len(seq_record.seq)
		
		print >> flankingcoords_file, orfID, beg, end
	predictfinal_file.close()
	flankingcoords_file.close()
	getoutput("tigr-glimmer extract -w '%s/%s.fna'  'tmpbia/%s.flankingcoords'    >   '%s/orfs_%s.flanking.fna' "  % (seq_record.id, seq_record.id, seq_record.id, seq_record.id, seq_record.id))	

	#Extraccion de regiones upstream de los nuevos codones start, para posterior busqueda de TFBS con MAST
	getoutput("upstream-coords.awk  500  0  'tmpbia/%sfinal.coords'  > 'tmpbia/orfs_%s.upstream.coords'   " % (file_prefix, file_prefix) )	
	getoutput("tigr-glimmer extract  -w  '%s/%s.fna' 'tmpbia/orfs_%s.upstream.coords'  > '%s/orfs_%s.upstream.fna'  "  % (seq_record.id, seq_record.id, file_prefix, seq_record.id, seq_record.id))	
f_contigs.close()

