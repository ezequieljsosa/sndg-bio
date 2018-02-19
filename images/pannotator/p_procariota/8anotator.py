#!/usr/bin/python

import sys
import os
import time
import getopt
from Bio import SeqIO
from commands import getoutput
import re

def help():
	print  "Anotacion final del genoma.\n\
Requiere los archivos de salida de los scripts previos\n\
Opciones:\n\
	-i	Archivo del genoma (o de proteoma, si solo se anota desde el paso 3). Default: contigs.fasta\n\
	-o	Nombre del organismo. Default: Bacteria\n\
	-l	Locus tag para los genes. Default: BIA\n\
	-s	Prioridad de anotacion complementaria (0: no usar, 1: usar si no hay otro dato, 2: usar si no hay datos Uniref, 3: usar si no hay datos Tigrfam, 4: primera opcion. Default: 0\n\
	-h 	Imprime este mensaje de ayuda\n\
        -p      datos de hmm. Default /data/pfamtigrfam/shorthmm	\n\
        -t      mapeo tigr-go. Default /data/pfamtigrfam/TIGRFAMS_GO_LINK	\n\
        -g      directorio de nombres de gen de uniprot. Default /data/uniprot/gsm/	\n\
        -m      directorio de mapeos de go-uniprot. Default /data/uniprot/goa	\n\
	-j      directorio de informacion de tigrfam. Default /data/pfamtigrfam/INFO \n\
	-v      mapeo pfam-go. Default /data/pfamtigrfam/pfam2go.txt \n\
	"


# Escapado de caracteres con significado especifico en gff3
def formatogff(cadena):
	cadenagff = re.sub("%", "%25", cadena, flags=re.I)
	cadenagff = re.sub(";", "%3B", cadenagff, flags=re.I)
	cadenagff = re.sub("=", "%3D", cadenagff, flags=re.I)
	cadenagff = re.sub("&", "%26", cadenagff, flags=re.I)
	cadenagff = re.sub(",", "%2C", cadenagff, flags=re.I)
	cadenagff = re.sub("\t", "%09", cadenagff, flags=re.I)
	cadenagff = re.sub("\n", "%0A", cadenagff, flags=re.I)
	cadenagff = re.sub("\r", "%0D", cadenagff, flags=re.I)
	return cadenagff

def sec(product, gene, symbol_source, product_source, EC, GO):
	if getoutput("grep %s  ./sec_%s/praze.out.btab" % (ORF.id, input_name)):
		if getoutput("grep %s  ./sec_%s/praze.out.btab" % (ORF.id, input_name)) and os.path.exists("sec_%s/praze.out.btab" % (input_name)) and float(getoutput("grep %s  ./sec_%s/praze.out.btab" % (ORF.id, input_name)).split("\n")[0].split("\t")[10]) >= 40:
			bstart, bend, bhstart, bhend = getoutput("grep %s  ./sec_%s/praze.out.btab" % (ORF.id, input_name)).split("\n")[0].split("\t")[6:10]
			if not paso_3:
				qstart, qend = ORF.description.strip().split()[2:4]
			else:
				qstart, qend = "1", str(len(ORF.seq) * 3)
			hstart, plus, hend, endline = getoutput("grep %s  ./sec_%s/praze.out.btab" % (ORF.id, input_name)).split("\n")[0].split("\t")[16:]
			cov = (float(bend) - float(bstart) + 1) / (abs(float(qend) - float(qstart)) + 1)    # cobertura en el query
			covh = (float(bhend) - float(bhstart) + 1) / (abs(float(hend) - float(hstart)) + 1)	 # cobertura en el hit
			cluster, uniprot = getoutput("grep %s  ./sec_%s/praze.out.btab" % (ORF.id, input_name)).split("\n")[0].split("\t")[5].split("_")
			evidence = 5

			if getoutput("grep %s -m 1 %s/genesymbolmap.%s" % (uniprot,params["g"] , uniprot[0:2])):
				gene, evidence = getoutput("grep %s -m 1 %s/genesymbolmap.%s" % (uniprot,params["g"] , uniprot[0:2])).split()[1:]
				symbol_source = cluster + "_" + uniprot
			if int(evidence) < 5:
				if cov >= 0.8 and covh >= 0.8:
					GOlines = getoutput("egrep \"%s\" '%/gp_association.goa_uniprot.%s'  " % (uniprot, params["m"] , uniprot[0:3])).split("\n")
					GO = list(GO)
					for line in GOlines:
						if line and line.split("\t")[3] not in GO:
							GO.append((line.split("\t")[3]))
					GO = ",".join(GO)
					product = getoutput("grep %s  ./sec_%s/praze.out.btab" % (ORF.id, input_name)).split("\n")[0].strip().split("\t")[15].split("n=")[0]	# + "     UNIREF100"
					product_source = cluster + "_" +  uniprot
				else:
					if covh < 0.8:
						GO = "GO:0008150,GO:0003674,GO:0005575"
					else:
						GOlines = getoutput("egrep \"%s\" '%s/gp_association.goa_uniprot.%s'  " % (uniprot,params["m"] , uniprot[0:3])).split("\n")
						GO = list(GO)
						for line in GOlines:
							if line and line.split("\t")[3] not in GO:
								GO.append((line.split("\t")[3]))
						GO = ",".join(GO)
					product = getoutput("grep %s  ./sec_%s/praze.out.btab" % (ORF.id, input_name) ).split("\n")[0].strip().split("\t")[15].split("n=")[0] + " domain protein"    # +"    UNIREF100"
					product_source = cluster + "_" + uniprot
					
	return product, gene, symbol_source, product_source, EC, GO

try:
	options, args = getopt.getopt(sys.argv[1:], "i:o:l:s:p:t:g:m:j:h")
except getopt.GetoptError as err:
	print str(err)
	sys.exit(2)


paso_3 = False

#Parametros configurables:
params = {}
params["i"] = "genoma.fasta"			# Archivo fasta del genoma
params["o"] = "Bacteria"				# Organismo
params["l"] = "BIA"						# Locus tag prefix
params["s"] = "0"						# Prioridad de datos complementarios
params["p"] = "/data/pfamtigrfam/shorthmm"			# Prioridad de datos complementarios
params["t"] = "/data/pfamtigrfam/TIGRFAMS_GO_LINK"
params["g"] = "/data/uniprot/gsm/"
params["m"] = "/data/uniprot/goa/"
params["j"] = "/data/pfamtigrfam/INFO"
params["v"] = "/data/pfamtigrfam/pfam2go.txt"
     
      
#Asigno los parametros que hayan sido definidos externamente
for option, value in options:
	if option.startswith("-"):
		option = option[1:]
	if option in params.keys():
		params[option] = value
	if option == "h":
		help()
		sys.exit()

log =  params["i"] + ".log"

sbt = "ncbi.sbt"

warning = ""
if not os.path.isfile("ncbi.sbt"):
	script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))
	getoutput("cp  %s/%s  ." %  (script_dir,sbt))

if not os.path.isfile(log):
	print "Advertencia. Esta iniciando el pipeline desde el Paso 8."
	warning = "True"
else:
	if int(getoutput("grep \"Paso 1\" %s -c" % log)) == 0 and int(getoutput("grep \"Paso 2\" %s -c" % log)) == 0 and int(getoutput("grep \"Paso 3\" %s -c" % log)) >= 1:
		paso_3 = True
		pass
	else:
		if int(getoutput("grep \"Paso 1\" %s -c" % log)) == 0:
			print "Advertencia: el Paso 1 no parece haber sido ejecutado."
			warning = "True"
		if int(getoutput("grep \"Paso 2\" %s -c" % log)) == 0:
			print "Advertencia: el Paso 2 no parece haber sido ejecutado."
			warning = "True"
		if int(getoutput("grep \"Paso 3\" %s -c" % log)) == 0:
			print "Advertencia: el Paso 3 no parece haber sido ejecutado."
			warning = "True"
		if int(getoutput("grep \"Paso 4\" %s -c" % log)) == 0:
			print "Advertencia: el Paso 4 no parece haber sido ejecutado."
			warning = "True"

	if int(getoutput("grep \"Paso 5\" %s -c" % log)) == 0:
		print "Advertencia: el Paso 5 no parece haber sido ejecutado."
		warning = "True"
	if int(getoutput("grep \"Paso 6\" %s -c" % log)) == 0:
		print "Advertencia: el Paso 6 no parece haber sido ejecutado."
		warning = "True"
	if int(getoutput("grep \"Paso 7\" %s -c" % log)) == 0:
		print "Advertencia: el Paso 7 no parece haber sido ejecutado."
		warning = "True"
#	if int(getoutput("grep \"Paso 8\" %s -c" % log)) >= 1:
#		print "Advertencia: ya se ejecuto el Paso 8.\n"
#		warning = "True"

if warning:
	print "Continuar de todos modos (s/n)?"
	respuesta = raw_input()
	if respuesta == "s":
		pass
	else:
		sys.exit(0)

log_step = open(log, "a")
if re.search("\.faa$", params["i"]):
	paso_3 = True

#Chequeo la existencia de los archivos
if getoutput("ls '%s' " % params["i"]) != params["i"]:
	print ("Error: No se encuentra el archivo %s \nUse el parametro -h para obtener ayuda sobre el uso del script" % params["i"])
	sys.exit(2)

if not paso_3:
	if not os.path.exists("ncbi.sbt"):
		print ("Error: No se encuentra el archivo ncbi.sbt")
		sys.exit(2)


rank = {"equivalog": 1, "ber": 2, "subfamily": 3, "equivalog_domain": 4, "superfamily": 5, "subfamily_domain": 6, "domain": 7, "pfam": 8, "tmhmm": 9,
        "lipoprotein": 10, "hypoth_equivalog": 11, "hypoth_equivalog_domain": 1000, "signature": 1000, "repeat": 1000, "paralog": 1000, "paralog_domain":1000, "exception": 1000}

SO_term = { "lipoprotein_signal_peptide":"SO:0100009", "signal_peptide":"SO:0000418", "transmembrane":"SO:0001077", "polypeptide_domain":"SO:0000417","beta_strand":"SO:0001111","peptide_helix":"SO:0001114", "disordered_region":"SO:0100003"}

structure = { "E":"beta_strand", "H":"peptide_helix","^":"disordered_region", "-":"disordered_region"}



if paso_3 == True:
	mock_contig = params["i"].split(".")[0];
	params["i"] = mock_contig + ".fna";

f_contigs = open(params["i"])
for seq_record in SeqIO.parse(f_contigs, "fasta"):
	print "Anotando", seq_record.id

	input_name = seq_record.id
	input_number = input_name[-3:]
	if paso_3:
		input_number = ""
		input_name = seq_record.id
	#Apertura y parseado de los  archivos de rRNAs, tRNAs, ORFs y salidas de Blast
	if not paso_3:
		rawgenome = open("%s/%s.fna" % (seq_record.id, seq_record.id))
		genome_iterator = SeqIO.parse(rawgenome, "fasta")
		rRNAs = open("%s/rRNA_%s.fna" % (seq_record.id, seq_record.id))
		rRNA_iterator = SeqIO.parse(rRNAs, "fasta")
		tRNAs = open("%s/tRNA_%s.fna" % (seq_record.id, seq_record.id))
		tRNA_iterator = SeqIO.parse(tRNAs, "fasta")
		TBL = open("%s/ncbi_%s.tbl" % (seq_record.id, seq_record.id), "w")
		orfs_fasta_final = open("%s/ncbi_%s.faa" % (seq_record.id, seq_record.id), "w")
		curatedhmmer_out = open("%s/curatedorfs_%s.hmmer" % (seq_record.id, seq_record.id))
		curatedhmmer_lines = curatedhmmer_out.readlines()
	
	else:
		toy_genome = "toy_genome.fasta"
		rawgenome = open(toy_genome, "w")
		print >> rawgenome, ">Chomosome_N"
		print >> rawgenome, ">ACTGACTGACTGACTG"
		rawgenome.close()
		rawgenome = open(toy_genome)
		genome_iterator = SeqIO.parse(rawgenome, "fasta")
		getoutput("mkdir %s/ber_%s " % (seq_record.id, input_name))
		
		getoutput("tail -n +2 %s/curatednr_%s.blast | gawk -F '	'  '{print $1\"	fecha	some_number	praze	path_source	\"$2\"	\"$4\"	\"$5*3\"	\"$6\"	\"$7\"	\"$8\"	similarity	false_cov	X	Y	\"$3\"	1	Plus	\"$5*100/$11\"	\"}'   > %s/ber_%s/praze.out.btab"  % (seq_record.id, input_name, seq_record.id, input_name))
		getoutput("mkdir %s/bercurated_%s " %  (seq_record.id, input_name))
		getoutput("tail -n +2 %s/curatednr_%s.blast | gawk -F '	'  '{print $1\"	fecha	some_number	praze	path_source	\"$2\"	\"$4\"	\"$5*3\"	\"$6\"	\"$7\"	\"$8\"	similarity	false_cov	X	Y	\"$3\"	1	Plus	\"$5*100/$11\"	\"}'   > %s/bercurated_%s/praze.out.btab"  % (seq_record.id, input_name, seq_record.id, input_name))

	if os.path.exists("%s/sec_%s.blast" % (seq_record.id, seq_record.id)):
		getoutput("mkdir %s/sec_%s " % (seq_record.id, seq_record.id))
		getoutput("tail -n +2 %s/sec_%s.blast | gawk -F '	'  '{print $1\"	fecha	some_number	praze	path_source	\"$2\"	\"$4\"	\"$5*3\"	\"$6\"	\"$7\"	\"$8\"	similarity	false_cov	X	Y	\"$3\"	1	Plus	\"$5*100/$11\"	\"}'   > %s/sec_%s/praze.out.btab"  % (seq_record.id, seq_record.id, seq_record.id, seq_record.id))

	ORFs = open("%s/orfs_%s.faa" % (seq_record.id, seq_record.id))
	ORF_iterator = SeqIO.parse(ORFs, "fasta")
	blast_out_cog = open("%s/bbhcog_%s.blast" % (seq_record.id, seq_record.id))
	blast_lines_cog = blast_out_cog.readlines()
	blast_out_ec = open("%s/ec_%s.blast" % (seq_record.id, seq_record.id))
	blast_lines_ec = blast_out_ec.readlines()
	try:
		psipred_out = open("%s/orfs_%s.psipred" % (seq_record.id, seq_record.id))	
		psipred_lines = psipred_out.readlines()
		disopred_out = open("%s/orfs_%s.disopred" % (seq_record.id, seq_record.id))	
		disopred_lines = disopred_out.readlines()
	except:
		psipred_lines = []
		disopred_lines = []
	# blast_out_int= open("interevidence_%s.blast" % input_name)
	# blast_lines_int = blast_out_int.readlines()
	hmmer_out = open("%s/orfs_%s.hmmer" % (seq_record.id, seq_record.id))
	hmmer_lines = hmmer_out.readlines()
		
	# psscan_out = open("orfs_%s.psscan" % input_name)
	# psscan_lines = psscan_out.readlines()
	tmhmm_out = open("%s/orfs_%s.tmhmm" % (seq_record.id, seq_record.id))
	tmhmm_lines = tmhmm_out.readlines()
	
	GFF3 = open("%s/ncbi_%s.gff3.features" % (seq_record.id, seq_record.id), "w")
	GFF3fasta = open("%s/ncbi_%s.gff3.fasta" % (seq_record.id, seq_record.id), "w")
	
	gff3id = "_".join(params["o"].split())
	print "\tPreprocesamiento..."	
	j = 0
	dict_blast_lines_cog = {}
	while j < len(blast_lines_cog):
		id = blast_lines_cog[j].split()[0]
		dict_blast_lines_cog[id] = blast_lines_cog[j].split("\t")[2].strip()
		print" %.2f" % float(100*j/len(blast_lines_cog))+ "% COG\r",
		j += 1
	
	j = 1
	best_blast_lines_ec = {}
	while j < len(blast_lines_ec):
		id = blast_lines_ec[j].split()[0]
		best_blast_lines_ec[id] = blast_lines_ec[j].split("\t")[2].strip()
		print" %.2f" % float(100*j/len(blast_lines_ec))+ "% EC \r",
		while j < len(blast_lines_ec) and id == blast_lines_ec[j].split()[0]:
			j += 1

	j = 0
	psipred = {}
	while j < len(psipred_lines):	
		if re.match("orf", psipred_lines[j]):

			psi_orf = psipred_lines[j].strip()
			psipred[psi_orf] = []
			j += 1
			base_preds = ""
			while j < len(psipred_lines) and not re.match("orf", psipred_lines[j]):		
				if re.match("Pred", psipred_lines[j]):
					base_preds = base_preds + psipred_lines[j].split()[1].strip()
					j += 1
				else:
					j += 1	
			for pred in re.finditer("H+|E+", base_preds):
				psipred[psi_orf].append([pred.start() + 1, pred.end(), pred.group(0)[0]])
		else:
			j += 1
			
	j = 0
	disopred = {}
	while j < len(disopred_lines):	
		if re.match("orf", disopred_lines[j]):
			diso_orf = disopred_lines[j].strip()
			disopred[diso_orf] = []
			j += 1
			base_preds = ""
			while j < len(disopred_lines) and not re.match("orf", disopred_lines[j]):		
				base_preds = base_preds + disopred_lines[j].strip().split()[2].strip()
				j += 1
			for pred in re.finditer("[\^\-]+", base_preds):
				disopred[diso_orf].append([pred.start() + 1, pred.end(), pred.group(0)[0]])
		else:
			j += 1
	
	
	# j=1
	# best_blast_lines_int =  []
	# while j < len(blast_lines_int):
		# id = blast_lines_int[j].split()[0]
		# best_blast_lines_int.append(blast_lines_int[j].strip())
		# while   j < len(blast_lines_int) and id == blast_lines_int[j].split()[0] :
			# j+=1

	j = 1
	tigr_class = {}
	pfam_match = {}
	if not paso_3:	
		while j < len(curatedhmmer_lines):
			print" %.2f" % float(100*j/len(curatedhmmer_lines))+ "% HMM refinados\r",	
			id = curatedhmmer_lines[j].split("\t")[0]
			isotype = 100
			gene_sym = ""
			EC_number = ""
			score_pfam = 0
			while   j < len(curatedhmmer_lines) and id == curatedhmmer_lines[j].split("\t")[0]:
				hitid = curatedhmmer_lines[j].split("\t")[2].strip()
				score = float(curatedhmmer_lines[j].split("\t")[-1])

				if score >= float(getoutput("egrep -A 14 '\s\s%s$'  %s" % (hitid,params["p"] )).split("\n")[13].strip().split()[-1][0:-1]):   # Si el score supera el trusted cutoff
					if re.match("TIGR", hitid):
						newisotype = rank[getoutput("egrep \"^IT\" '%s/%s.INFO'  " % (params["j"],hitid)).split()[-1]]
					else:
						newisotype = rank["pfam"]
						covq = float(curatedhmmer_lines[j].split("\t")[-3].strip())
						covh = 100 * abs(float(curatedhmmer_lines[j].split("\t")[-5].strip()) - float(curatedhmmer_lines[j].split("\t")[-6].strip())) / float(getoutput("egrep -A 4 '\s\s%s$'  %s" % (hitid,params["p"] )).split("\n")[3].strip().split()[-1])
						# print id, covq, covh
						if getoutput("egrep -A 7 '\s\s%s$'  %s" % (hitid,params["p"] )).split("\n")[6].strip().split()[-1] == "no":
							if covq < 30 or covh < 85:
								newisotype = 100
		
						beg, end, hb, he, ident, cov, ev, score = curatedhmmer_lines[j].split("\t")[4:]
						if float(score) > score_pfam:
							score_pfam = float(score)
							pfam_id = getoutput("egrep -A 1 '\s\s%s$'  %s" % (hitid,params["p"] )).split("\n")[1].strip().split()[1]
							pfam_match[id] = [ curatedhmmer_lines[j].split("\t")[2].strip(), score_pfam, beg, end, pfam_id]      # nombre, hitid, beg, end, score
					if newisotype < isotype:
						isotype = newisotype
						if isotype != 8:
							gene_sym = getoutput("egrep \"^GS\" '%s/%s.INFO'  " % (params["j"],hitid))
							if gene_sym:
								gene_sym = gene_sym.split()[-1]
							EC_number = getoutput("egrep \"^EC\" '%s/%s.INFO'  " % (params["j"],hitid))
							if EC_number:
								EC_number = EC_number.split()[-1]
						tigr_class[id] = curatedhmmer_lines[j].split("\t")[3].strip(), gene_sym, EC_number, isotype, hitid   # nombre, gene_sym, EC_number, isotype, hitid
		
				j += 1
	
	j = 1
	
	while j < len(hmmer_lines):
		print" %.2f" % float(100*j/len(hmmer_lines))+ "% HMM restantes \r",		
		id = hmmer_lines[j].split("\t")[0]
		if not tigr_class.has_key(id):
			isotype = 100
			gene_sym = ""
			EC_number = ""
			score_pfam = 0
			pfam_haskey = "True"
			while j < len(hmmer_lines) and id == hmmer_lines[j].split("\t")[0]:
				hitid = hmmer_lines[j].split("\t")[2].strip()
				score = float(hmmer_lines[j].split("\t")[-1])
				if score >= float(getoutput("egrep -A 14 '\s\s%s$'  '%s'  " % (hitid,params["p"] )).split("\n")[13].strip().split()[-1][0:-1]):
					if re.match("TIGR", hitid):
						newisotype = rank[getoutput("egrep \"^IT\" '%s/%s.INFO'  " % (params["j"],hitid)).split()[-1]]
					else:
						newisotype = rank["pfam"]
						covq = float(hmmer_lines[j].split("\t")[-3].strip())
						covh = 100*abs(float(hmmer_lines[j].split("\t")[-5].strip()) - float(hmmer_lines[j].split("\t")[-6].strip()))/float(getoutput("egrep -A 4 '\s\s%s$'  %s" % (hitid,params["p"])   ).split("\n")[3].strip().split()[-1])
						if   getoutput("egrep -A 7 '\s\s%s$'  %s" % (hitid,params["p"] )   ).split("\n")[6].strip().split()[-1] == "no":
							if  covq < 30 or covh < 85:
							# newisotype = 100
								newisotype = 100
	
						beg, end, hb, he, ident, cov, ev, score = hmmer_lines[j].split("\t")[4:]
						if float(score) > score_pfam:
							score_pfam = float(score)
							pfam_id = getoutput("egrep -A 1 '\s\s%s$'  %s" % (hitid,params["p"] )).split("\n")[1].strip().split()[1]
							if not pfam_match.has_key(id):
								pfam_haskey = "False"	
							if pfam_haskey == "False":
								pfam_match[id] = [hmmer_lines[j].split("\t")[2].strip(), score_pfam, beg, end, pfam_id]      # nombre, hitid, beg, end, score
					if newisotype <  isotype:
						isotype = newisotype
						if isotype !=8:
							gene_sym =  getoutput("egrep \"^GS\" '%s/%s.INFO'  " % (params["j"],hitid) )
							if gene_sym:
								gene_sym = gene_sym.split()[-1]
							EC_number = getoutput("egrep \"^EC\" '%s/%s.INFO'  " % (params["j"],hitid ))
							if EC_number:
								EC_number = EC_number.split()[-1] 					
						tigr_class[ id ] = hmmer_lines[j].split("\t")[3].strip(), gene_sym, EC_number, isotype, hitid       #nombre, gene_sym, EC_number, isotype, hitid
				j += 1
		else:
			j += 1

	#Formateado del archivo fasta del genoma
	genome = genome_iterator.next()

	if not paso_3:
		print >> TBL, ">Features " + genome.description

	print >> GFF3,  "##gff-version 3\n##feature-ontology so.obo"

	if not paso_3:
		gffcontig =  formatogff(input_name) + "\tv0.1-bia_pap\tcontig\t1\t" + str(len(genome.seq)) + "\t.\t+\t.\tID=" + formatogff(input_name) + ";translation_table=11;organism_name=" + formatogff(gff3id) + ";abbreviation=" + formatogff(params["l"]) + ";defline=" + formatogff(genome.description) + ";Name=" + formatogff(input_name)
		print >> GFF3, gffcontig
		genome.description = genome.id + " [organism=" + params["o"] + "] " + genome.description[len(genome.id)+1:]
		new_genome = open("%s/ncbi_%s.fna" % (seq_record.id, seq_record.id), "w")
		SeqIO.write(genome, new_genome, "fasta")
		genome.description = ""
		SeqIO.write(genome, GFF3fasta, "fasta")

	#Generacion del archivo .sbt

	#Generacion de archivo .tbl
	#ORFs predichos
	print "\tProteinas..."
	i = 1
	len_cds = int(getoutput("grep '>' -c %s/orfs_%s.faa" % (seq_record.id, seq_record.id)))
	while True:
		#print" %.2f" % float(100*i/len_cds)+ "% CDS\r",
		product = ""
		gene = ""
		symbol_source = ""
		product_source = ""
		EC = ""
		GO = ""
		try:
			ORF = ORF_iterator.next()
		except StopIteration:
			break
		if params["s"] == "4":
			product, gene, symbol_source, product_source, EC, GO = sec(product, gene, symbol_source, product_source, EC, GO)

		if (product == "" or re.search("hypothetical protein", product)) and tigr_class.has_key(ORF.id):
			if tigr_class[ORF.id][3] == 1:
				gene = tigr_class[ORF.id][1]
				if gene:
					symbol_source = tigr_class[ORF.id][4]
				EC = tigr_class[ORF.id][2]
				GOlines = getoutput("egrep \"%s\" '%s'  " % (tigr_class[ORF.id][4], params["t"])).split("\n")
				for line in GOlines:
					if line:
						GO = GO + line.split("\t")[1] + ","
				product = tigr_class[ORF.id][0].split(": ")[1]    # + "    TIGRFAM"
				product_source = tigr_class[ORF.id][4]
	
		if params["s"] == "3":
			if (product == "" or re.search("hypothetical protein", product)):
				product, gene, symbol_source, product_source, EC, GO = sec(product, gene, symbol_source, product_source, EC, GO)
		
		if product == "" and getoutput("grep %s  %s/ber_%s/praze.out.btab" % (ORF.id, seq_record.id,seq_record.id)):
			if getoutput("grep %s  %s/bercurated_%s/praze.out.btab" % (ORF.id, seq_record.id, seq_record.id)) and os.path.exists("%s/bercurated_%s/praze.out.btab" % (seq_record.id, seq_record.id)) and float(getoutput("grep %s  %s/bercurated_%s/praze.out.btab" % (ORF.id, seq_record.id, seq_record.id)).split("\n")[0].split("\t")[10]) >= 40:
				bstart, bend, bhstart, bhend = getoutput("grep %s  %s/bercurated_%s/praze.out.btab" % (ORF.id, seq_record.id, seq_record.id)).split("\n")[0].split("\t")[6:10]
				if not paso_3:
					qstart, qend = ORF.description.strip().split()[2:4]
				else:
					qstart, qend = "1", str(len(ORF.seq) * 3)
				hstart, plus, hend, endline = getoutput("grep %s  %s/bercurated_%s/praze.out.btab" % (ORF.id, seq_record.id, seq_record.id)).split("\n")[0].split("\t")[16:]
				cov = (float(bend) - float(bstart) + 1) / (abs(float(qend) - float(qstart)) + 1)    # cobertura en el query
				covh = (float(bhend) - float(bhstart) + 1) / (abs(float(hend) - float(hstart)) + 1)	 # cobertura en el hit
				cluster, uniprot = getoutput("grep %s  %s/bercurated_%s/praze.out.btab" % (ORF.id, seq_record.id, seq_record.id)).split("\n")[0].split("\t")[5].split("_")
				evidence = 5
				if getoutput("grep %s -m 1 %s/genesymbolmap.%s" % (uniprot,params["g"] , uniprot[0:2])):
					gene, evidence = getoutput("grep %s -m 1 %s/genesymbolmap.%s" % (uniprot,params["g"] , uniprot[0:2])).split()[1:]
					symbol_source = cluster + "_" + uniprot
				if int(evidence) < 5:
					if cov >= 0.8 and covh >= 0.8:
						GOlines = getoutput("egrep \"%s\" '%s/gp_association.goa_uniprot.%s'  " % (uniprot,params["m"], uniprot[0:3])).split("\n")
						GO = list(GO)
						for line in GOlines:
							if line and line.split("\t")[3] not in GO:
								GO.append((line.split("\t")[3]))
						GO = ",".join(GO)
						product = getoutput("grep %s  %s/bercurated_%s/praze.out.btab" % (ORF.id, seq_record.id, seq_record.id)).split("\n")[0].strip().split("\t")[15].split("n=")[0]	# + "     UNIREF100"
						product_source = cluster + "_" + uniprot
					else:
						if covh < 0.8:
							GO = "GO:0008150,GO:0003674,GO:0005575"
						else:
							GOlines = getoutput("egrep \"%s\" '%s/gp_association.goa_uniprot.%s'  " % (uniprot,params["m"], uniprot[0:3])).split("\n")
							GO = list(GO)
							for line in GOlines:
								if line and line.split("\t")[3] not in GO:
									GO.append((line.split("\t")[3]))
							GO = ",".join(GO)
						product = getoutput("grep %s  %s/bercurated_%s/praze.out.btab" % (ORF.id, seq_record.id, seq_record.id) ).split("\n")[0].strip().split("\t")[15].split("n=")[0] + " domain protein"    # +"    UNIREF100"
						product_source = cluster + "_" + uniprot
			
			elif getoutput("grep %s  %s/ber_%s/praze.out.btab" % (ORF.id, seq_record.id, seq_record.id) ) and os.path.exists("%s/ber_%s/praze.out.btab" % (seq_record.id, seq_record.id)) and float(getoutput("grep %s  %s/ber_%s/praze.out.btab" % (ORF.id, seq_record.id, seq_record.id) ).split("\n")[0].split("\t")[10]) >= 40:
				if float(getoutput("grep %s  %s/ber_%s/praze.out.btab" % (ORF.id, seq_record.id, seq_record.id)).split("\n")[0].split("\t")[10]) >= 40:
					bstart, bend, bhstart, bhend = getoutput("grep %s  %s/ber_%s/praze.out.btab" % (ORF.id, seq_record.id, seq_record.id)).split("\n")[0].split("\t")[6:10]
					qstart, qend = ORF.description.strip().split()[2:4]
					hstart, plus, hend, endline = getoutput("grep %s  %s/ber_%s/praze.out.btab" % (ORF.id, seq_record.id, seq_record.id)).split("\n")[0].split("\t")[16:]
					cov = (float(bend) - float(bstart) + 1) / (abs(float(qend) - float(qstart)) + 1)
					covh = (float(bhend) - float(bhstart) + 1) / (abs(float(hend) - float(hstart)) + 1)	 # cobertura en el hit
					cluster, uniprot = getoutput("grep %s  %s/ber_%s/praze.out.btab" % (ORF.id, seq_record.id, seq_record.id)).split("\n")[0].split("\t")[5].split("_")
					evidence = 5
					if getoutput("grep %s -m 1 %s/genesymbolmap.%s" % (uniprot,params["g"], uniprot[0:2])):
						gene, evidence = getoutput("grep %s -m 1 %s/genesymbolmap.%s" % (uniprot,params["g"], uniprot[0:2])).split()[1:]
						symbol_source = cluster + "_" + uniprot
					if int(evidence) < 5:
						if cov >= 0.8 and covh >= 0.8:
							#print getoutput("egrep \"%s\" '/data/uniprot/goa/gp_association.goa_uniprot.%s'  " % (uniprot, uniprot[0:3])).split("\n")
							GOlines = getoutput("egrep \"%s\" '%s/gp_association.goa_uniprot.%s'  " % (uniprot,params["m"], uniprot[0:3])).split("\n")
							GO = list(GO)
							for line in GOlines:
								if line and line.split("\t")[3] not in GO:
									GO.append(  (line.split("\t")[3]) )
							GO = ",".join(GO)
							product = getoutput("grep %s  %s/ber_%s/praze.out.btab" % (ORF.id, seq_record.id, seq_record.id) ).split("\n")[0].strip().split("\t")[15].split("n=")[0]              #+  "     UNIREF100"
							product_source = cluster + "_" + uniprot
						else:
							if covh < 0.8:
								GO = "GO:0008150,GO:0003674,GO:0005575"
							else:
								GOlines =  getoutput("egrep \"%s\" '%s/gp_association.goa_uniprot.%s'  " % (uniprot, params["m"], uniprot[0:3])   ).split("\n")
								GO= list(GO)
								for line in GOlines:
									if line and line.split("\t")[3] not in GO:
										GO.append(  (line.split("\t")[3]) )
								GO = ",".join(GO)				
							product =   getoutput("grep %s  %s/ber_%s/praze.out.btab" % (ORF.id, seq_record.id, seq_record.id) ).split("\n")[0].strip().split("\t")[15].split("n=")[0]+ " domain protein"   #+"     UNIREF100"

							product_source = cluster + "_" + uniprot

		if params["s"] == "2":
			if (product == "" or re.search("hypothetical protein", product)):
				product, gene, symbol_source, product_source, EC, GO = sec(product, gene, symbol_source, product_source, EC, GO)
	
		if  product == ""  and tigr_class.has_key(ORF.id):
			EC = tigr_class[ ORF.id][2]
			if tigr_class[ ORF.id ][3] == 3:
				GOlines =  getoutput("egrep \"%s\" '%s'  " % (tigr_class[ ORF.id][4],params["t"])).split("\n")
				for line in GOlines:
					if line:
						GO = GO + line.split("\t")[1] + ","
				product = tigr_class[ ORF.id][0].split(": ")[1] + " family protein"  # +"    TIGRFAM"
				product_source = tigr_class[ ORF.id][4]
			elif tigr_class[ ORF.id ][3] == 4:
				GOlines =  getoutput("egrep \"%s\" '%s'  " % (tigr_class[ ORF.id][4],params["t"])).split("\n")
				for line in GOlines:
					if line:
						GO = GO + line.split("\t")[1] + ","
				product = tigr_class[ ORF.id][0].split(": ")[1]  +"    TIGRFAM"
				product_source = tigr_class[ ORF.id][4]
			elif tigr_class[ ORF.id ][3] == 5:
				GOlines =  getoutput("egrep \"%s\" '%s'  " % (tigr_class[ ORF.id][4],params["t"])).split("\n")
				for line in GOlines:
					if line:
						GO = GO + line.split("\t")[1] + ","
				product = tigr_class[ ORF.id][0].split(": ")[1] + " family protein"      #+"    TIGRFAM"
				product_source = tigr_class[ ORF.id][4]
			elif tigr_class[ ORF.id ][3] == 6:
				GOlines =  getoutput("egrep \"%s\" '%s'  " % (tigr_class[ ORF.id][4],params["t"])).split("\n")
				for line in GOlines:
					if line:
						GO = GO + line.split("\t")[1] + ","
				product = tigr_class[ ORF.id][0].split(": ")[1] + " domain protein"     #+"    TIGRFAM"
				product_source = tigr_class[ ORF.id][4]
			elif tigr_class[ ORF.id ][3] == 7:
				GOlines =  getoutput("egrep \"%s\" '%s'  " % (tigr_class[ ORF.id][4],params["t"])).split("\n")
				for line in GOlines:
					if line:
						GO = GO + line.split("\t")[1] + ","
				product = tigr_class[ ORF.id][0].split(": ")[1] + " domain protein"       #+"    TIGRFAM"
				product_source = tigr_class[ ORF.id][4]
			elif tigr_class[ ORF.id ][3] == 8:
				GOlines =  getoutput("egrep \" %s \" '%s'  " % (tigr_class[ ORF.id][4], params["v"]) ).split("\n")
				for line in GOlines:
					if line:
						GO = GO + line.split()[-1] + ","
				product = tigr_class[ORF.id][0] + " family protein"      #+"    PFAM" 	
				product_source = tigr_class[ ORF.id][4]
			elif tigr_class[ ORF.id ][3] == 11:
				GOlines =  getoutput("egrep \"%s\" '%s'  " % (tigr_class[ ORF.id][4],params["t"]) ).split("\n")
				for line in GOlines:
					if line:
						GO = GO + line.split("\t")[1] + ","
				product = "conserved hypothetical protein"			
				product_source = tigr_class[ ORF.id][4]
		tmhmmpred = getoutput("grep %s  %s/orfs_%s.tmhmm" % (ORF.id, seq_record.id, seq_record.id))
			
		if product == "" or (tigr_class.has_key(ORF.id) and tigr_class[ORF.id][3] == 11):
			if tmhmmpred and float(tmhmmpred.split("\t")[2].split("=")[1]) /  float(tmhmmpred.split("\t")[1].split("=")[1]) >= 0.5:
	
				GO = "GO:0008150,GO:0003674,GO:0016020"
				product = "putative membrane protein"  #+"TMHMM"
				product_source = "TMHMM"
		if product == "" or (tigr_class.has_key(ORF.id) and tigr_class[ ORF.id ][3] == 11):
			lipopred = getoutput("grep %s  %s/orfs_%s.lipop" % (ORF.id, seq_record.id, seq_record.id))
			if lipopred:
				lipopred = lipopred.split("\n")[0].split()[2].strip()
				lipomargin = float(getoutput("grep %s  %s/orfs_%s.lipop" % (ORF.id, seq_record.id, seq_record.id)).split("\n")[0].split()[4].split("=")[1].strip())
			if lipopred == "SpII" and lipomargin >= 3:
				GO = "GO:0008150,GO:0003674,GO:0016020"
				product = "putative lipoprotein" #+"    LIPOP"
				product_source = "lipoP"

			signalpred = getoutput("grep %s  %s/orfs_%s.signalp" % (ORF.id, seq_record.id, seq_record.id))
			if signalpred:
				product = "putative signal peptide protein"
				product_source = "signalP"
		if params["s"] == "1":
			if (product == "" or re.search("hypothetical protein", product)):
				product, gene, symbol_source, product_source, EC, GO = sec(product, gene, symbol_source, product_source, EC, GO)
	
		if product == "":
			GO = "GO:0008150,GO:0003674,GO:0005575"
			product = "hypothetical protein"   #+ "     DEFAULT"
			product_source = "default"
		product=re.sub("superfamily family", "family", product, flags=re.I)
		product=re.sub("family[\W]*family", "family", product, flags=re.I)
		product=re.sub("domain[\W]*domain", "domain", product, flags=re.I)
		product=re.sub("domain family protein", "domain protein", product, flags=re.I)
		if re.search("putative", product, flags=re.I) and  re.search("family", product, flags=re.I):
			product=re.sub("putative[\W]", "", product, flags=re.I)
		product=re.sub("^the[\W]", "", product, flags=re.I)
		if re.search("PAS domain|S-box domain protein", product, flags=re.I):
			product= "sensory box protein"
		if re.search("response regulator", product, flags=re.I):
			product= "response regulator"
		if re.search("conserved hypothetical|protein.*unknown function", product, flags=re.I):
			product= "conserved hypothetical protein"
		if re.search("DUF.*protein|ncharacterized ACR", product):
			product= "conserved hypothetical protein"
		if re.search("Y\S*[A-Z] family", product):
			product= "conserved hypothetical protein"
		product=re.sub("subunit family protein", "subunit", product, flags=re.I)
		if re.search("protein \S+ protein", product, flags=re.I):
			product = re.sub("protein ", "", product, count=1, flags=re.I)
		if re.search("family \S+ family", product, flags=re.I):
			product = re.sub("family", "", product, count=1, flags=re.I)
		product=re.sub(" or ", "/", product, flags=re.I)
		product=re.sub("^predicted\s+", "", product, flags=re.I)
		product=re.sub("precursor\s+", "", product, flags=re.I)
		product=re.sub("probable", "putative", product, flags=re.I)
		product=re.sub("^\s*", "", product, flags=re.I)
		product=re.sub("\(ec.+?\)", "", product, flags=re.I)
		if re.search("^[A-Z][a-z]", product):
			product= product[0].lower()+ product[1:]
		product=re.sub("\s+", " ", product, flags=re.I)
		product=re.sub("\s+$", "", product, flags=re.I)
		product=re.sub("\.$", "", product, flags=re.I)  	
		
	   	if not paso_3:
			beg, end = ORF.description.strip().split()[2:4]
		else:
		    beg, end = "1", str(len(ORF.seq) * 3)
		if gene == "":
			if getoutput("grep %s  %s/ber_%s/praze.out.btab" % (ORF.id, seq_record.id, seq_record.id)) and os.path.exists("%s/ber_%s/praze.out.btab" % (seq_record.id, seq_record.id)):
				cluster, uniprotid = getoutput("grep %s  %s/ber_%s/praze.out.btab" % (ORF.id, seq_record.id, seq_record.id)).split("\n")[0].split("\t")[5].split("_")
				if getoutput("grep %s -m 1 %s/genesymbolmap.%s" % (uniprotid,params["g"], uniprotid[0:2]) ):	
					gene =   getoutput("grep %s -m 1 %s/genesymbolmap.%s" % (uniprotid,params["g"], uniprotid[0:2]) ).split()[1]
					symbol_source =  cluster + "_" + uniprotid
		if not paso_3:
			print >> TBL, "%s	%s	gene" % (beg, end)
			print >> TBL, "			locus_tag	%s_%s" % (params["l"], i)
			if gene:
				print >> TBL, "			gene	%s"  % gene
			print >> TBL, "%s	%s	CDS" % (beg, end)
			print >> TBL, "			protein_id	gnl|bia|%s_%s" % (params["l"], i)
	
		if re.search("domain|family",product) and not re.search("protein",product):
			product = product + " protein"
		if not paso_3:
			print >> TBL, "			product	%s" % product
			# print >> TBL, "			transcript_id"
			if dict_blast_lines_cog.has_key(ORF.id) and dict_blast_lines_cog[ORF.id]:
				print >> TBL, "			note	%s" % dict_blast_lines_cog[ORF.id]
			if EC:
				print >> TBL, "			EC_number	%s" % EC
			elif best_blast_lines_ec.has_key(ORF.id):
				print >> TBL, "			EC_number	%s" % best_blast_lines_ec[ORF.id]	
	
		if int(beg) < int(end):
			dir = "+"
		else:
			dir = "-"
			beg, end = end, beg
			
		db_ref = ""
		if paso_3:
			header_mapping_file = params["i"] + ".header_mapping"
			if os.path.exists(header_mapping_file):
				db_ref = ";Alias=" + getoutput("grep %s %s" % (ORF.id, header_mapping_file)).split("\t")[0].strip()
			else:
				db_ref = ";Alias=" + ORF.id
			
		gffgene =   formatogff(input_name) + "\tv0.1bia_pap\tgene\t" + beg + "\t" + end + "\t.\t"+ dir + "\t.\tID=" + formatogff(gff3id) + ".gene." + input_number + str(i) + db_ref
		gffexon =   formatogff(input_name) + "\tv0.1bia_pap\texon\t" + beg + "\t" + end + "\t.\t"+dir+"\t0\tID=" + formatogff(gff3id) + ".exon." + input_number +str(i) + ";Parent=" + formatogff(gff3id )+ ".transcript." + input_number +str(i)
		gffmRNA =   formatogff(input_name) + "\tv0.1bia_pap\tmRNA\t" + beg + "\t" + end + "\t.\t"+dir+"\t.\tID=" + formatogff(gff3id) + ".transcript." + input_number + str(i) + ";Parent=" + formatogff(gff3id) + ".gene." + input_number + str(i)
		if symbol_source:
			gffmRNA = gffmRNA + ";gene_symbol_source=" + symbol_source
		if GO:
			gffmRNA = gffmRNA + ";GO=" + GO 	
		gffmRNA = gffmRNA + ";gene_product_name_source=" + product_source
		if EC:
			gffmRNA = gffmRNA + ";EC=" + re.sub(";", ",", EC, flags=re.I)
		elif best_blast_lines_ec.has_key(ORF.id):
			gffmRNA = gffmRNA + ";EC=" + re.sub(";", ",", best_blast_lines_ec[ORF.id],  flags=re.I)
		if gene:
			gffmRNA = gffmRNA + ";gene_symbol=" + gene
	#	if gene2:
	#		print "Gene2 ",
	#		if gene != gene2:		
	#			print "distinto"
	#			gffmRNA = gffmRNA + "," + gene2
	
		gffmRNA = gffmRNA + ";Note=" + product + ";description=" + product +  ";Name=" + formatogff(gff3id) + ".transcript." + input_number + str(i) + ";locus_tag=" + params["l"] + "_" +  str(i)
		gffCDS =   formatogff(input_name) + "\tv0.1bia_pap\tCDS\t" + beg + "\t" + end + "\t.\t"+dir+"\t0\tID=" + formatogff(gff3id) + ".CDS." + input_number + str(i) + ";Parent=" + formatogff(gff3id) + ".transcript." + input_number + str(i)
		if dict_blast_lines_cog.has_key(ORF.id) and dict_blast_lines_cog[ORF.id]:
			gffmRNA = gffmRNA + ";top_cog_hit=" + formatogff(dict_blast_lines_cog[ORF.id])
	
		gffpfam = ""
		if pfam_match.has_key(ORF.id):
			if dir == '+':
				dom_beg = int(pfam_match[ORF.id][2]) + int(beg) - 1
				dom_end = int(pfam_match[ORF.id][3]) + int(beg) - 1
			if dir == '-':
				dom_end = int(end) - int(pfam_match[ORF.id][2]) + 1
				dom_beg = int(end) - int(pfam_match[ORF.id][3]) + 1
			
	
			gffpfam = formatogff(input_name) + "\tv0.1bia_pap\tpolypeptide_domain\t" + str(dom_beg) + "\t" + str(dom_end) + "\t" + str(pfam_match[ORF.id][1]) + "\t"+ dir +"\t." + "\tID=" + formatogff(gff3id) + ".domain." + input_number + str(i) + ";Parent=" + formatogff(gff3id) + ".gene." + input_number + str(i) + ";Ontology_term=" + formatogff(SO_term["polypeptide_domain"]) +  ";description=" + str(pfam_match[ORF.id][0]) + ";Accession=" + str(pfam_match[ORF.id][4])
	
		print >> GFF3, gffgene
		print >> GFF3, gffexon
		print >> GFF3, gffmRNA
		print >> GFF3, gffCDS
		print >> GFF3, gffpfam
	
		gfftmhmm = ""
		if tmhmmpred:
			topology = tmhmmpred.split("\t")[5].split("=")[1]
			topology_clean = topology.replace("o", "i")
			top_array = topology_clean.split("i")[1:-1]
			h = 0
			for segment in top_array:
				h = h + 1
				seg_begp, seg_endp = segment.split("-")
				if dir == '+':
					seg_beg = int(seg_begp) + int(beg) - 1
					seg_end = int(seg_endp) + int(beg) - 1
				if dir == '-':
					seg_beg = int(end) - int(seg_endp) + 1
					seg_end = int(end) - int(seg_begp) + 1
					
	
				gfftmhmm = formatogff(input_name) + "\tv0.1bia_pap\ttransmembrane\t" + str(seg_beg) + "\t" + str(seg_end) + "\t.\t"+ dir +"\t." + "\tID=" + formatogff(gff3id) + ".transmembrane." + input_number + str(i) + "." + str(h) + ";Parent=" + formatogff(gff3id) + ".gene." + input_number + str(i)  + ";Ontology_term=" + formatogff(SO_term["transmembrane"])
	
				print >> GFF3, gfftmhmm
				
	
		gfflipop = ""
		lipopred = getoutput("grep %s  %s/orfs_%s.lipop" % (ORF.id, seq_record.id, seq_record.id))		
		if lipopred:
			lipopredT = lipopred.split("\n")[0].split()[2].strip()
			lipomargin = float(getoutput("grep %s  %s/orfs_%s.lipop" % (ORF.id, seq_record.id, seq_record.id)).split("\n")[0].split()[4].split("=")[1].strip())	
			if lipopredT == "SpII" and lipomargin >= 3:
				cleavbegp, cleavendp = lipopred.split("\n")[0].split()[5].split("=")[1].split("-")[:]
				if dir == '+':
					cleavbeg = int(beg)
					cleavend = int(cleavbegp) + int(beg) - 1
				if dir == '-':
					cleavend = int(end)
					cleavbeg = int(end) - int(cleavbegp) + 1
				gfflipop = formatogff(input_name) + "\tv0.1bia_pap\tlipoprotein_signal_peptide\t" + str(cleavbeg) + "\t" + str(cleavend) + "\t.\t"+ dir +"\t." + "\tID=" + formatogff(gff3id) + ".lipoprotein." + input_number + str(i) + ";Parent=" + formatogff(gff3id) + ".gene." + input_number + str(i) + ";Ontology_term=" + formatogff(SO_term["lipoprotein_signal_peptide"])
	
				print >> GFF3, gfflipop				
	
		gffsignalp = ""
		signalpred = getoutput("grep %s  %s/orfs_%s.signalp" % (ORF.id, seq_record.id, seq_record.id))
		if signalpred:
			cleavbegp = getoutput("grep %s  %s/orfs_%s.signalp" % (ORF.id, seq_record.id, seq_record.id)).split()[6]
			if dir == '+':
				cleavbeg = int(beg)
				cleavend = int(cleavbegp) + int(beg) - 1
			if dir == '-':
				cleavend = int(end)			
				cleavbeg = int(end) - int(cleavbegp) + 1
	
			gffsignalp = formatogff(input_name) + "\tv0.1bia_pap\tsignal_peptide\t" + str(cleavbeg) + "\t" + str(cleavend) + "\t.\t"+ dir +"\t." + "\tID=" + formatogff(gff3id) + ".signalpeptide." + input_number + str(i) + ";Parent=" + formatogff(gff3id) + ".gene." + input_number + str(i) + ";Ontology_term=" + formatogff(SO_term["signal_peptide"])
			print >> GFF3, gffsignalp	
	
				
		if ORF.id in psipred:
			print "feo"
			h = 0
			for region in psipred[ORF.id]:
				h = h + 1
				seg_begp, seg_endp, stype = region[:]
				if dir == '+':
					seg_beg = int(seg_begp) + int(beg) - 1
					seg_end = int(seg_endp) + int(beg) - 1
				if dir == '-':
					seg_beg = int(end) - int(seg_endp) + 1
					seg_end = int(end) - int(seg_begp) + 1
					
	
				gffpsipred = formatogff(input_name) + "\tv0.1bia_pap\t" + structure[stype] + "\t" + str(seg_beg) + "\t" + str(seg_end) + "\t.\t"+ dir +"\t." + "\tID=" + formatogff(gff3id) + "." + structure[stype] + "." + input_number + str(i) + "." + str(h) + ";Parent=" + formatogff(gff3id) + ".gene." + input_number + str(i)  + ";Ontology_term=" + formatogff(SO_term[structure[stype]])
				print >> GFF3, gffpsipred 
		
		if ORF.id in disopred:
			print "feo"
			h = 0
			for region in disopred[ORF.id]:
				h = h + 1
				seg_begp, seg_endp, stype = region[:]
				if dir == '+':
					seg_beg = int(seg_begp) + int(beg) - 1
					seg_end = int(seg_endp) + int(beg) - 1
				if dir == '-':
					seg_beg = int(end) - int(seg_endp) + 1
					seg_end = int(end) - int(seg_begp) + 1
					
			
				gffdisopred = formatogff(input_name) + "\tv0.1bia_pap\t" + structure[stype] + "\t" + str(seg_beg) + "\t" + str(seg_end) + "\t.\t"+ dir +"\t." + "\tID=" + formatogff(gff3id) + "." + structure[stype] + "." + input_number + str(i) + "." + str(h) + ";Parent=" + formatogff(gff3id) + ".gene." + input_number + str(i)  + ";Ontology_term=" + formatogff(SO_term[structure[stype]])
				print SO_term[structure[stype]]
				print structure[stype]
				print >> GFF3, gffdisopred 
							
	
		if not paso_3:
			ORF.id = formatogff(gff3id) + ".CDS." + input_number + str(i)
			ORF.description = ""
			SeqIO.write(ORF, GFF3fasta, "fasta")
			ORF.id = "bia|" + ORF.id + "_" + genome.id + " " + product + " OS=" + params["o"] + " GN=" + gene + " GO=" + GO
			SeqIO.write(ORF, orfs_fasta_final, "fasta")
	
		i += 1
	
	
	# #ORFs encontrados en regiones de interevidencia
	# j=0
	# while True:
		# try:
			# ORFdesc, IDhit, product, dfrom, dto, tfrom, tto, ident, cov, evalue  = best_blast_lines_int[j].split("\t")
			# ORFid = ORFdesc.split()[0]
			# if re.match("hypothetical",product):
				# product = "hypothetical protein"
			# j+=1
		# except:
			# break
	
		# print >> TBL, "%s	%s	gene" % (dfrom, dto)
		# print >> TBL, "			locus_tag	BIA_%s" % i
		# print >> TBL, "%s	%s	CDS" % (dfrom, dto)
		# print >> TBL, "			protein_id	gnl|bia|BIA_%s" % i		
		# if re.search("domain|family",product) and not re.search("protein",product):
			# product = product + " protein"
		# print >> TBL, "			product	%s" % product
		# # print >> TBL, "			transcript_id"
		# if dict_blast_lines_cog.has_key(ORFid):
			# print >> TBL, "			note	%s" %dict_blast_lines_cog[ORFid]
		# if best_blast_lines_ec.has_key(ORF.id):
			# print >> TBL, "			EC_number	%s" % best_blast_lines_ec[ORFid]	
		# i+=1	
	

	if not paso_3:
		len_trna = int(getoutput("grep '>' -c %s/tRNA_%s.fna"  %  (seq_record.id, seq_record.id) ))
		len_rrna = int(getoutput("grep '>' -c %s/rRNA_%s.fna"  % (seq_record.id, seq_record.id) ))
		#tRNAs predichos
		print "\ttRNAs..."
		while True:
			try:
		    		tRNA = tRNA_iterator.next()	
				print" %.2f" % float(100*i/len_trna)+ "% tRNA\r",				
		    	except StopIteration:
				break
	    	
		    	beg, end = tRNA.id.strip().split("_")[-1].split("-")
		    	tRNA_product = tRNA.description.split("/molecule=")[1].split()[0]
		    	anticodon =  tRNA.description.split("/anticodon=")[1]
		    	tRNAscore =  tRNA.description.split("/score=")[1].split()[0]
		    	print >> TBL, "%s	%s	gene" % (beg, end)
		    	print >> TBL, "			locus_tag	%s_%s" % (params["l"], i)
		    	print >> TBL, "%s	%s	tRNA" % (beg, end)
		    	print >> TBL, "			product	%s" % tRNA_product
	    	
		    	if int(beg) < int(end):
		    		dir = "+"
		    	else:
		    		dir = "-"
		    		beg, end = end, beg
		    	gffgene =   formatogff(input_name) + "\tv0.1bia_pap\tgene\t" + beg + "\t" + end + "\t.\t"+dir+"\t.\tID=" + formatogff(gff3id) + ".gene." + input_number + str(i)
		    	gffexon =   formatogff(input_name) + "\tv0.1bia_pap\texon\t" + beg + "\t" + end + "\t.\t"+dir+"\t0\tID=" + formatogff(gff3id) + ".exon." +  input_number + str(i)
		    	gfftRNA =   formatogff(input_name) + "\tv0.1bia_pap\ttRNA\t" + beg + "\t" + end + "\t.\t"+dir+"\t.\tID=" + formatogff(gff3id) + ".tRNA." +  input_number + str(i) + ";tRNA_anti-codon=" + anticodon + ";score=" + tRNAscore +  ";Note=" + tRNA_product + ";description=" + tRNA_product + ";locus_tag=" + params["l"] + "_" +  str(i)
		    	gffCDS =   formatogff(input_name) + "\tv0.1bia_pap\tCDS\t" + beg + "\t" + end + "\t.\t"+dir+"\t0\tID=" + formatogff(gff3id) + ".CDS." +  input_number + str(i)
		    	print >> GFF3, gffgene
		    	print >> GFF3, gffexon
		    	print >> GFF3, gfftRNA
		    	print >> GFF3, gffCDS	
		    	i+=1
	
		#rRNAs predichos
		print "\trRNAs..."
		while True:
			try:
				rRNA = rRNA_iterator.next()	
				print" %.2f" % float(100*i/len_rrna)+ "% rRNA\r",	
			except StopIteration:
				break
			
			begend, dir = rRNA.id.strip().split("_")[-2:]
			beg, end = begend.split("-")
			if dir == "DIR-":
				beg, end = end, beg
				dir = "-"
			else:
				dir = "+"
			rRNA_product = rRNA.description.split("/molecule=")[1].split()[0]
			rRNAscore =  rRNA.description.split("/score=")[1]
			print >> TBL, "%s	%s	gene" % (beg, end)
			print >> TBL, "			locus_tag	%s_%s" % (params["l"], i)
			print >> TBL, "%s	%s	rRNA" % (beg, end)
			print >> TBL, "			product	%s" % rRNA_product	
	
			gffgene = formatogff(input_name) + "\tv0.1bia_pap\tgene\t" + beg + "\t" + end + "\t.\t"+dir+"\t.\tID=" + formatogff(gff3id) + ".gene." + input_number + str(i)
			gffexon = formatogff(input_name) + "\tv0.1bia_pap\texon\t" + beg + "\t" + end + "\t.\t"+dir+"\t0\tID=" + formatogff(gff3id) + ".exon." +  input_number + str(i)
			gfftRNA = formatogff(input_name) + "\tv0.1bia_pap\trRNA\t" + beg + "\t" + end + "\t.\t"+dir+"\t.\tID=" + formatogff(gff3id) + ".rRNA." +  input_number + str(i) +  ";score=" + rRNAscore +  ";Note=" + rRNA_product + ";description=" + rRNA_product + ";locus_tag=" + params["l"] + "_" +  str(i)
			gffCDS =   formatogff(input_name) + "\tv0.1bia_pap\tCDS\t" + beg + "\t" + end + "\t.\t"+dir+"\t0\tID=" + formatogff(gff3id) + ".CDS." +  input_number + str(i)
			print >> GFF3, gffgene
			print >> GFF3, gffexon
			print >> GFF3, gfftRNA
			print >> GFF3, gffCDS	
			
			i += 1
	

	print >> GFF3, "##FASTA"

	if not paso_3:
		rRNAs.close()
		tRNAs.close()
		TBL.close()
		new_genome.close()
		orfs_fasta_final.close()
		curatedhmmer_out.close()
	rawgenome.close()
	ORFs.close()
	blast_out_cog.close()
	blast_out_ec.close()
	# blast_out_int.close()
	hmmer_out.close()
	# psscan_out.close()
	tmhmm_out.close()
	
	GFF3.close()
	GFF3fasta.close()
	


	#Generacion de archivo .gbf
	print "\tGenerando archivos..."
	if not paso_3:	
		print getoutput("tbl2asn -j \"[gcode=11]\" -k m  -Z discrepancies  -t ncbi.sbt -V vb -i %s/ncbi_%s.fna  -f  %s/ncbi_%s.tbl  -k m" % (seq_record.id, seq_record.id, seq_record.id, seq_record.id))


	#Generacion de archivo .gff3
	print getoutput("cat %s/ncbi_%s.gff3.features > %s/ncbi_%s.gff3" % (seq_record.id, seq_record.id, seq_record.id, seq_record.id) )
	if not paso_3:
		print getoutput("cat %s/ncbi_%s.gff3.fasta >> %s/ncbi_%s.gff3" % (seq_record.id, seq_record.id, seq_record.id, seq_record.id) )
	getoutput("rm  %s/ncbi_%s.gff3.features && rm %s/ncbi_%s.gff3.fasta" % (seq_record.id, seq_record.id, seq_record.id, seq_record.id) )
	
now = time.strftime("%d/%m/%Y %H:%M:%S", time.gmtime())
print >> log_step, "[ %s ] Paso 8 Anotator corrido sobre archivo %s. Organismo: %s. Locus tag para los genes: %s" % (now, params["i"], params["o"], params["l"] )
print >> log_step, "Pipeline finalizado."

getoutput("rm toy_genome.fasta")
log_step.close()
