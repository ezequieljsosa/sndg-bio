#!/usr/bin/python


#   genbank_to_tbl.py "my organism name" "my strain ID" "ncbi project id" < my_sequence.gbk
#   writes seq.fsa, seq.tbl as output

import sys
import re
import getopt
from copy import copy
from Bio import SeqIO

def help():
    print  "Conversion tbl -> gff3.\n\
Opciones:\n\
	-i	Archivo de entrada tbf. Default: contigs.tbl\n\
	-t	Archivo de salida gff3. Default: archivo_de_entrada.tbl\n\
	-f	Archivo de salida fasta. Default: archivo_de_entrada.fasta\n\
	-o	Nombre del organismo Default: vacio\n\
	-s	ID de Cepa. Default: vacio\n\
	-n	ID de NCBI Project. Default: vacio\n\
	-h 	Imprime este mensaje de ayuda\n"

try:
    options, args = getopt.getopt(sys.argv[1:], "i:t:f:o:s:n:h")
except getopt.GetoptError as err:
    print str(err)
    sys.exit(2)

params = {}
params["i"] = "contigs.tbl"                   	
params["t"] = ""
params["f"] = ""
params["o"] = ""
params["s"] = ""
params["n"] = ""
params["e"] = "1e-5"                   			

for option, value in options:
    if option.startswith("-"): option = option[1:]
    if option in params.keys(): params[option] = value
    if option == "h":
        help()
        sys.exit()

if not params["t"]:
    params["t"] = ".".join(params["i"].split(".")[:-1]) + ".tbl"
if not params["f"]:
    params ["f"] = ".".join(params["i"].split(".")[:-1]) + ".fasta"

def formatogff(cadena):
	cadenagff=re.sub("%", "%25", cadena, flags=re.I)
	cadenagff=re.sub(";", "%3B", cadenagff, flags=re.I)
	cadenagff=re.sub("=", "%3D", cadenagff, flags=re.I)
	cadenagff=re.sub("&", "%26", cadenagff, flags=re.I)
	cadenagff=re.sub(",", "%2C", cadenagff, flags=re.I)
	cadenagff=re.sub("\t", "%09", cadenagff, flags=re.I)
	cadenagff=re.sub("\n", "%0A", cadenagff, flags=re.I)
	cadenagff=re.sub("\r", "%0D", cadenagff, flags=re.I)
	return cadenagff


def find_gene_entry(features, locus_tag):
    for f in features:
        if f.type == 'gene':
            if f.qualifiers['locus_tag'][0] == locus_tag:
                return f
    print locus_tag
    raise ValueError

coding = ['CDS', 'tRNA', 'rRNA', 'ncRNA']

seqid = 0
fasta_fh = open("seq.fsa", "w")
feature_fh = open("seq.tbl", "w")
allowed_tags = ['locus_tag', 'gene', 'product', 'pseudo', 'protein_id', 'gene_desc', 'old_locus_tag']
records = list(SeqIO.parse(sys.stdin, "genbank"))

for rec in records:
    for f in rec.features:
        if f.type in coding and 'gene' in f.qualifiers:
            print f.qualifiers['locus_tag'][0]

            f2 = find_gene_entry(rec.features, f.qualifiers['locus_tag'][0])
            f2.qualifiers['gene'] = f.qualifiers['gene']

            del f.qualifiers['gene']

for rec in records:
    seqid += 1

    if len(rec) <= 200:
        print >>sys.stderr, "skipping small contig %s" % (rec.id,)
        continue

    circular = rec.annotations.get('molecule', 'linear')
    rec.description = "[organism=%s] [strain=%s] [topology=%s] [molecule=DNA] [tech=wgs] [gcode=11]" % (sys.argv[1], sys.argv[2], circular)
    SeqIO.write([rec], fasta_fh, "fasta")

    print >>feature_fh, ">Feature %s" % (rec.name,)
    for f in rec.features:
        if f.strand == 1:
            print >>feature_fh, "%d\t%d\t%s" % (f.location.nofuzzy_start + 1, f.location.nofuzzy_end, f.type)
        else:
            print >>feature_fh, "%d\t%d\t%s" % (f.location.nofuzzy_end, f.location.nofuzzy_start + 1, f.type)

        if f.type == 'CDS' and 'product' not in f.qualifiers:
            f.qualifiers['product'] = ['hypothetical protein']

        if f.type == 'CDS':
            f.qualifiers['protein_id'] = ["gnl|ProjectID_%s|%s" % (sys.argv[3], f.qualifiers['locus_tag'][0])]

        if f.type in coding:
            del f.qualifiers['locus_tag']

        for key, vals in f.qualifiers.iteritems():
            my_allowed_tags = copy(allowed_tags)
            if 'pseudo' in f.qualifiers:
                my_allowed_tags.append('note')

            if key not in my_allowed_tags:
                continue

            for v in vals:
                if len(v) or key == 'pseudo':
                    print >>feature_fh, "\t\t\t%s\t%s" % (key, v)

