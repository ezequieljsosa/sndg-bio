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
	-d 	Base de datos general para busqueda de similitud. Default: data/uniprot/uniref100/uniref100.fasta\n\
	-s	Subset de Base de datos Uniprot de organismos similares al target. Default: vacio\n\
	-e	Umbral de e-value para la busqueda de similitud. Default: 1e-5\n\
	-h 	Imprime este mensaje de ayuda\n"

try:
    options, args = getopt.getopt(sys.argv[1:], "i:d:e:s:h")
except getopt.GetoptError as err:
    print str(err)
    sys.exit(2)

#Parametros configurables:	
params = {}
params["i"] =  "contigs.fasta"                   	        #Archivo de contigs de entrada por defecto
params["d"] = "/data/uniprot/uniref100/uniref100.fasta" 	#Base de datos general por defecto
params["s"] = ""                                 		#Base de datos por defecto
params["e"] = "1e-5"                   				#Umbral de e-value

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

# Uso de Blast y BER para Busqueda de posibles errores de secuenciacion que pudieran malograr la prediccion de genes
# (debido a frameshifts o falsos codon de stop generados por sustitucion)
# Los genes que superan el chequeo mantienen los hits de blast para la primera anotacion.
print "-Anotacion preliminar con Blast"
getoutput("mkdir tmpbia")

if params["s"]:
    print "    Busqueda en organismos relacionados (%s)\n" % params["s"].split("/")[-1]
    base_inicial = params["s"]

f_contigs = open(params["i"])
for seq_record in SeqIO.parse(f_contigs, "fasta"):
	print "	Procesando",  seq_record.id
	#Parseo del archivo fasta con los ORFs predichos inicialmente
	queryfile = open("%s/orfs_%s.faa" % (seq_record.id, seq_record.id) )
	queryrecs = list(SeqIO.parse(queryfile, "fasta"))
	queryfile.close()
	# Ejecucion de blastp

	getoutput("blastall -p blastp -i '%s/orfs_%s.faa' -d %s -e %s -a 4 -v 20 -b 20 -o 'tmpbia/nr_%s.blast' " % (seq_record.id, seq_record.id,base_inicial,params["e"],seq_record.id))

	#Parseo del archivo blast generado
	fileoutHandler = open("tmpbia/nr_"+ seq_record.id+".blast")
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
				cov = int(100*(float(hsp.query_end)-float(hsp.query_start)+1.0)/len(queryrec.seq))  #Calculo la cobertura
				ident = int(float(hsp.identities[0])/float(hsp.identities[1])*100)    #Porcentaje de identidad
				score = float(hsp.score) 	#Score
				evalue = hsp.expect          	#E-value
				dfrom = int(hsp.query_start)	#Posicion de inicio y fin del hsp en el query
				dto = int(hsp.query_end)
				tfrom = int(hsp.sbjct_start)  		#Posicion de inicio y fin del hsp en el hit
				tto = int(hsp.sbjct_end)
				#Genero una lista con los datos del parseo para cada alineamiento
				hits.append((queryrec.id, queryrec.description,sbjctID,sbjttProd,dfrom,dto,tfrom,tto,ident,cov,evalue,score))

	fileoutHandler.close()

	#Salida del ID del query, ID y cadena del pdb, posiciones de los hsp, % de identidad, cobertura y e-value
	blast_out_nr = open("%s/nr_%s.blast" % (seq_record.id, seq_record.id), "w")
	print >> blast_out_nr, "ID del query        Descripcion        ID del hit        Producto        Hsp del query        Hsp del hit        % de identidad        Cobertura        e-value"
	for hit in hits:
		print >> blast_out_nr,  "\t".join([ str(f) for f in hit[1:-1]])
	blast_out_nr.close()
	
	print "	Refinamiento de genes predichos"
	#Ejecucion de BER
	startsites = {}	
	if len(hits)!=0:
		getoutput("formatdb -i   '%s/orfs_%s.flanking.fna'  -o T -p F  "   %  (seq_record.id, seq_record.id) )
		getoutput("ber -i 'tmpbia/nr_%s.blast' -d %s -f formatdb -D  '%s/orfs_%s.flanking.fna'   -F formatdb  -m  '%s/orfs_%s.map'     -e 1e-5 -E 1e-5  -n 20  -N 0"    %   (seq_record.id, base_inicial, seq_record.id, seq_record.id,seq_record.id, seq_record.id))
	
		#Analisis de overlaps y curado de sitios de start y stop.
		hmmev = open("%s/orfs_%s.hmmer" % (seq_record.id, seq_record.id))
		hmmev_lines = hmmev.readlines()
		hmmhits =  []
		for ev in hmmev_lines[1:]:
			hmmhits.append(ev.split("\t")[0])
		hmmev.close()
		blasthits = [ hit[0] for hit in hits ]
		overlapcutoff = 60
		identitycutoff = 60
		pvaluecutoff  = 1e-30
		characterizedbonusvote = 4
		minvotecutoff = 2	
		flank_startuse =  open("flank_startuse")
		flank, startuse = flank_startuse.readlines()
		flank_startuse.close()
		flank = int(flank)
		startuseatg, startusegtg, startusettg = startuse.strip().split(",")
		new_ber_lines =  []
		votos = 0
		for hit in hits:
			if float(hit[10]) <= pvaluecutoff:
				pot_new_ber_line = getoutput("grep  '%s'  out/praze.out.btab  |  grep '%s' " % (hit[0], hit[2])  )
				if  pot_new_ber_line and float(pot_new_ber_line.split("\t")[10] )  >= identitycutoff:
					new_ber_lines.append(pot_new_ber_line )		
				
		# for hit in hits:
			# pot_new_ber_line = getoutput("grep  '%s'  out/praze.out.btab  |  grep '%s' " % (hit[0], hit[2])  )
			# if float(hit[10]) <= pvaluecutoff or float(pot_new_ber_line.split("\t")[10] )  >= identitycutoff:
				# new_ber_lines.append( pot_new_ber_line )		
	
		#Parseo del archivo fasta con flanqueo de 300bp
		flankfile = open( "%s/orfs_%s.flanking.fna"   %  (seq_record.id, seq_record.id))
		flankrecs = SeqIO.to_dict(SeqIO.parse(flankfile, "fasta"))
		flankfile.close()
	
		#Busqueda de start sites en el mismo frame del hit
		print "	Curado de sitios de start"
		frame =  flank % 3

		for orf in (flankrecs.keys()):
			starts = []
			iterator = re.finditer("atg|gtg|ttg", str(flankrecs[orf].seq), flags=re.IGNORECASE)
			for match in iterator:
				if match.start() % 3 == frame:
					starts.append([match.start()+1,0] )
			startsites[orf] = starts
	
		#Cada start site recibe un voto si coincide con el start del hit de blast
		for ber_line in new_ber_lines:
			start = int(ber_line.split("\t")[6])
			orf =  ber_line.split("|")[1].split("\t")[0]
			for starts in startsites[orf][:]:
				if start == starts[0]:
					starts[1]+=1
					# if is_characterized: +4
		
		# Se retienen el o los start con valor maximo. Si no se llega al cutoff de votos, se elimina
		for orf in (flankrecs.keys()):
			startsites[orf].sort(key=lambda start: start[1], reverse=True)
			for starts in startsites[orf][:]:
				if starts[1] < minvotecutoff or starts[1] < startsites[orf][0][1] :
					startsites[orf].remove(starts)
			#Si hay empate, se verifican rbs, startuse, start original y secuencia mas larga, en ese orden 			
			if not startsites[orf]:
				del startsites[orf]	
			elif len(startsites[orf]) > 1:
				for starts in startsites[orf][:]:
					rbsrange = str(flankrecs[orf].seq)[starts[0]-21:starts[0]-6]
					for pos in range(0, 10):
						if float(len(re.findall("A|a|G|g", rbsrange[pos:pos+6]) ) )/6  >= 0.75:
							seqstart = str(flankrecs[orf].seq).lower()[starts[0]-1:starts[0]+2]
							if seqstart == "atg": starts[1]+=1+float(startuseatg)
							elif seqstart == "gtg": starts[1]+=1+float(startusegtg)
							elif seqstart == "ttg": starts[1]+=1+float(startusettg)
							break
				startsites[orf].sort(key=lambda start: start[1], reverse=True)
				for starts in startsites[orf][:]:
					if starts[1] < startsites[orf][0][1] :
						startsites[orf].remove(starts)
				if len(startsites[orf]) > 1:
					if flank+1  in  [start[0] for start in startsites[orf] ]:
						del startsites[orf]
					else:
						startsites[orf].sort()
						startsites[orf] = [startsites[orf][0]]
							
	# Extraccion final de ORFs verificando overlaps
	print "	Extraccion final de ORFs"
	rrnafile = open("%s/rRNA_%s.fna" % (seq_record.id, seq_record.id))
	trnafile = open("%s/tRNA_%s.fna" % (seq_record.id, seq_record.id))
	orffile = open("%s/orfs_%s.faa" % (seq_record.id, seq_record.id))	
	rRNAs = list(SeqIO.parse(rrnafile, "fasta"))
	tRNAs = list(SeqIO.parse(trnafile, "fasta"))
	ORFs = list(SeqIO.parse(orffile, "fasta"))
	rrnafile.close()
	trnafile.close()
	orffile.close()
	neworffile = open("%s/orfs_%s.faa" % (seq_record.id, seq_record.id), "w")	
	RNAs = []
	for rRNA in rRNAs:
		RNAs.append( (rRNA.id).split("_")[-2].split("-") )
	for tRNA in tRNAs:
		beg, end =  (tRNA.id).split("_")[-1].split("-")
		if int(beg) > int(end):
			beg, end = end, beg
		RNAs.append( [beg, end] )
	ORFsov = []
	for ORF in ORFs:
		orfid, contigID, beg, end, score =  ORF.description.strip().split()[:]
		if int(beg) > int(end):
			beg, end = end, beg
		ORFsov.append( [orfid, contigID, beg, end] )

	for ORF in ORFs:
		curated = ""
		orfid, contigID, beg, end, score =  ORF.description.strip().split()[:]
		beg, end = int(beg), int(end)
		if  startsites.has_key(orfid):
			if beg < end:
				beg =  beg + startsites[orfid][0][0] - (flank +1)
			else:
				beg =  beg - startsites[orfid][0][0] + (flank +1)				
			dnaseq =  Seq(str(flankrecs[orfid].seq)[startsites[orfid][0][0]-1:-flank], generic_dna)
			protseq = str(dnaseq.translate(table=11))
			curated = "curated"
		else:	
			protseq =  str(ORF.seq)
		overlap = 0
		if beg < end:
			for RNA in RNAs:
				if  beg<=int(RNA[1]) and end>=int(RNA[0]):
					overlap = 1
					break
			for ORFov in ORFsov:
				if  beg<=int(ORFov[3])- 60 and end>=int(ORFov[2])+60:
					if ORFov[0] in blasthits[:] or ORFov[0]  in hmmhits[:]:	
						overlap = 1
						break
		else:
			for RNA in RNAs:
				if  end<=int(RNA[1]) and beg>=int(RNA[0]):
					overlap = 1
					break	
			for ORFov in ORFsov:
				if  end<=int(ORFov[3])- 60 and beg>=int(ORFov[2]) + 60:
					if ORFov[0] in blasthits[:] or ORFov[0]  in hmmhits[:]:	
						overlap = 1
						break							
		if orfid in blasthits[:] or orfid in hmmhits[:] or not overlap: 				
			neworf=SeqRecord(Seq(protseq), id = orfid, description="%s %d %d %s %s" % (contigID, beg, end, score, curated))	
			SeqIO.write(neworf, neworffile, "fasta")
	neworffile.close()

	getoutput("mv out %s/ber_%s" % (seq_record.id, seq_record.id))
	
	#Enmascarado del genoma, para extraccion de regiones de interevidencia
	print "	Enmascarado del genoma"
	orffile = open("%s/orfs_%s.faa" % (seq_record.id, seq_record.id))	
	rrnafile = open("%s/rRNA_%s.fna" % (seq_record.id, seq_record.id))
	trnafile = open("%s/tRNA_%s.fna" % (seq_record.id, seq_record.id))
	rRNAs = list(SeqIO.parse(rrnafile, "fasta"))
	tRNAs = list(SeqIO.parse(trnafile, "fasta"))
	ORFs = list(SeqIO.parse(orffile, "fasta"))
	interevidence = open("%s/interevidence_%s.fna" % (seq_record.id, seq_record.id), "w")
	seq = list(seq_record.seq)
	contigID = seq_record.id

	#Enmascarado de genes
	for ORF in ORFs:
		beg, end=  ORF.description.strip().split()[2:4]
		beg, end = int(beg), int(end)
		if beg > end:
			beg, end = end, beg
		# beg, end = beg+50, end-50
		seq[beg-1:end] = ["N"]*(end-beg+1)
		seqname =  int(ORF.description.strip().split()[0][3:])

	#Enmascarado de rRNA
	for rRNA in rRNAs:
		beg, end=  rRNA.id.strip().split("_")[-2].split("-")
		beg, end = int(beg), int(end)
		# beg, end = beg+50, end-50
		seq[beg-1:end] = ["N"]*(end-beg+1)

	#Enmascarado de tRNA
	for tRNA in tRNAs:
		beg, end=  tRNA.id.strip().split("_")[-1].split("-")
		beg, end = int(beg), int(end)
		if beg > end:
			beg, end = end, beg
		# beg, end = beg+50, end-50
		seq[beg-1:end] = ["N"]*(end-beg+1)
	
	seq = "".join(seq)
	interev = re.split("[nN]{25,}", seq)
	pos = []
	for n in re.finditer("[nN]{25,}", seq):
		pos.append((n.start(), n.end()))	
	if len(pos)!=0 and pos[0][0] != 0:
		pos.insert(0,(0,0))
	else:
		del interev[0]
	j = 0
	while j<len(interev):
		beg = pos[j][1] + 1
		end = beg + len(interev[j])-1
		orf=SeqRecord(Seq(interev[j]), id = "orf%.5d" % (seqname+1+j), description="%s %d %d" % (contigID, beg, end))	
		if  len(interev[j]) >= 100:
			SeqIO.write(orf, interevidence, "fasta")
		j+=1
	
	rrnafile.close()
	trnafile.close()
	orffile.close()
	interevidence.close()	

	print "	Busqueda con Blastx en regiones de interevidencia\n\n"
	# Parseo del archivo de genoma con los ORFs enmascarados
	queryfile = open("%s/interevidence_%s.fna" % (seq_record.id, seq_record.id) )
	queryrecs = list(SeqIO.parse(queryfile, "fasta"))
	queryfile.close()
	## Ejecucion de blastx
	#getoutput("blastall -p blastx -i '%s/interevidence_%s.fna' -d %s -e %s -a 4 -v 20 -b 20 -o tmpbia/interevidence_%s.blast" % (seq_record.id, seq_record.id, base_inicial,params["e"], seq_record.id))
	
	#Parseo del archivo blast generado
	#fileoutHandler = open("tmpbia/interevidence_"+seq_record.id+".blast")
	#blastparser = BlastParser()
	#blast_iterator = Iterator(fileoutHandler, blastparser)

	#hits = []
	#i=0
	#for blast_record in blast_iterator:
	#	queryrec=queryrecs[i]
	#	i+=1
	#	for alignment in blast_record.alignments: 	#Para cada alineamiento de blastout
	#		sbjctID = alignment.title.split(" ")[0][1:]	#Obtengo ID de la secuencia subject
	#		sbjttProd = " ".join(alignment.title.split("[")[0].strip().split(" ")[1:])
	#		for hsp in alignment.hsps:                     		#Para cada hsp en cada alineamiento
	#			cov = abs(int(100*(float(hsp.query_end)-float(hsp.query_start)+1.0)/len(queryrec.seq)))  #Calculo la cobertura
	#			ident = int(float(hsp.identities[0])/float(hsp.identities[1])*100)    #Porcentaje de identidad
	#			score = float(hsp.score) 	#Score
	#			evalue = hsp.expect          	#E-value
	#			orfstart, orfend = queryrec.description.split()[2:]
	#			if int(orfstart) > int(orfend):
	#				orfstart = orfend
	#			dfrom = int(hsp.query_start) + int(orfstart) - 1	#Posicion de inicio y fin del hsp en el query
	#			dto = int(hsp.query_end) + int(orfstart) - 1
	#			tfrom = int(hsp.sbjct_start)  						#Posicion de inicio y fin del hsp en el hit
	#			tto = int(hsp.sbjct_end)
	#			#Genero una lista con los datos del parseo para cada alineamiento
	#			hits.append((queryrec.id, queryrec.description,sbjctID,sbjttProd,dfrom,dto,tfrom,tto,ident,cov,evalue,score))
	#fileoutHandler.close()
		
	#Salida del ID del query, ID y cadena del pdb, posiciones de los hsp, % de identidad, cobertura y e-value
	blast_out_int= open("%s/interevidence_%s.blast" % (seq_record.id, seq_record.id), "w")
	print >> blast_out_int, "ID del query        Descripcion        ID del hit        Producto        Hsp del query        Hsp del hit        % de identidad        Cobertura        e-value"
	#for hit in hits:
	#	print >> blast_out_int, "\t".join([ str(f) for f in hit[1:-1]])
	#blast_out_int.close()
			
getoutput("rm "+ "-r tmpbia")
f_contigs.close()

