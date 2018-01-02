#!/usr/bin/python


import sys
import re
import getopt
from copy import copy
from Bio import SeqIO


def help():
    print  "Conversion gbf -> tbl + gff3 + fasta.\n\
Opciones:\n\
	-i	Archivo de entrada gbf. Default: contigs.gbf\n\
	-t	Archivo de salida tbl. Default: archivo_de_entrada.tbl\n\
	-g	Archivo de salida gff3. Default: archivo_de_entrada.gff3\n\
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
params["i"] =  "contigs.gbf"                   	
params["t"] = ""
params["g"] = ""
params["f"] = ""
params["o"] = ""
params["s"] = ""
params["n"] = ""                  			

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
if not params["g"]:
    params ["g"] = ".".join(params["i"].split(".")[:-1]) + ".gff3"
if params["n"]:
    params["n"] = "_" + params["n"]

def find_gene_entry(features, locus_tag):
    for f in features:
        if f.type == 'gene':
            if f.qualifiers['locus_tag'][0] == locus_tag:
                return f
    raise ValueError

# Escapado de caracteres con significado especifico en gff3
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

coding = ['CDS', 'tRNA', 'rRNA']


seqid = 0
featid = 0
fasta_fh = open(params["f"], "w")
feature_fh = open(params["t"], "w")
gff3_fh = open(params["g"],"w")

allowed_tags = ['locus_tag', 'gene', 'product', 'pseudo', 'protein_id', 'gene_desc', 'old_locus_tag']
records = list(SeqIO.parse(params["i"], "genbank"))
gff3_featsCDS = ''



for rec in records:
    input_number = rec.name[-3:]
    for f in rec.features:
        if f.type in coding and 'gene' in f.qualifiers:
            f2 = find_gene_entry(rec.features, f.qualifiers['locus_tag'][0])
            f2.qualifiers['gene'] = f.qualifiers['gene']
            del f.qualifiers['gene']
    if 'locus_tag' in f.qualifiers:
        rec.locus_tag = f.qualifiers['locus_tag'][0].split("_")[0]
    else:
        rec.locus_tag = ""
    if 'transl_table' in f.qualifiers:
        rec.trans_tab = f.qualifiers['transl_table'][0].split("_")[0]
    else:
        rec.trans_tab = 11
    rec.defline = rec.name + " " + rec.description.split(".")[0]

for rec in records:
    seqid += 1
    mol_type = rec.annotations.get('molecule', 'circular')
    rec.description = "[organism=%s] [strain=%s] [topology=%s] [molecule=DNA] [tech=wgs] [gcode=11]" % (params["o"], params["s"], mol_type)
    SeqIO.write([rec], fasta_fh, "fasta")

    print >>feature_fh, ">Feature %s" % (rec.name)
    print >> gff3_fh, "##gff-version 3"
    print >> gff3_fh, "##feature-ontology so.obo"


    for f in rec.features:
        gff3_feats = {}
        if f.type == 'source':
            organism = f.qualifiers["organism"][0]
            gff3id = "_".join(organism.split())
        if f.strand == 1:
            start = f.location.nofuzzy_start + 1
            end = f.location.nofuzzy_end
            gff3_feats['strand'] = '+'
        else:
            start = f.location.nofuzzy_end
            end = f.location.nofuzzy_start + 1
            gff3_feats['strand'] = '-'

        print >>feature_fh, "%d\t%d\t%s" % (start, end, f.type)
        gff3_feats['start'] = f.location.nofuzzy_start + 1
        gff3_feats['end'] = f.location.nofuzzy_end

        if f.type == 'CDS' and 'product' not in f.qualifiers:
            f.qualifiers['product'] = ['hypothetical protein']

        if f.type == 'CDS':
            f.qualifiers['protein_id'] = ["gnl|ProjectID%s|%s" % (params["n"], f.qualifiers['locus_tag'][0])]

        if f.type == 'rRNA':

            f.qualifiers['product'] = [f.qualifiers['product'][0].split("S")[0] + "s_rRNA" ]

        if f.type in coding:
            del f.qualifiers['locus_tag']

        for key, vals in f.qualifiers.iteritems():
            my_allowed_tags = copy(allowed_tags)
            if 'pseudo' or 'note' in f.qualifiers:
                my_allowed_tags.append('note')
            if 'EC_number' in key:
                my_allowed_tags.append('EC_number')
                vals = [";".join(vals)]

            if key not in my_allowed_tags:
                continue

         #   print vals
            for v in vals:
                if len(v) or key == 'pseudo':
                    print >>feature_fh, "\t\t\t%s\t%s" % (key, v)
                    if key == 'gene':
                        gff3_feats['gene'] = v
                    if key == 'product':
                        gff3_feats['product'] = v
                    if key == 'EC_number':
                        gff3_feats['EC'] = ";EC=" + ",".join(v.split(";"))
                    if key == 'note':
                        if re.search("COG", v):
                            gff3_feats['note'] = ";top_cog_hit=" + formatogff(v)
                        if re.match("tRNA", v):
                            gff3_feats['note'] = v

        if f.type == 'source':
            print >> gff3_fh, "%s\tgenbankfile\tcontig\t%d\t%d\t.\t%s\t.\tID=%s;translation_table=%s;organism_name=%s;abbreviation=%s;defline=%s;Name=%s" % (formatogff(rec.id), gff3_feats['start'], gff3_feats['end'], gff3_feats['strand'], formatogff(rec.name), rec.trans_tab, formatogff(f.qualifiers["organism"][0]), formatogff(rec.locus_tag), formatogff(rec.defline),formatogff(rec.name))

        if  f.type == 'gene':
            if 'gene' in gff3_feats:
                gff3_featsCDS = ";gene_symbol=" + gff3_feats['gene']

            featid += 1
            print >> gff3_fh, "%s\tgenbankfile\tgene\t%d\t%d\t.\t%s\t.\tID=%s.gene.%s%s" % (formatogff(rec.id), gff3_feats['start'], gff3_feats['end'], gff3_feats['strand'], formatogff(gff3id), input_number, featid)

        if  f.type == 'CDS':
            if not 'EC' in gff3_feats:
                gff3_feats['EC'] = ''
            if not 'product' in gff3_feats:
                gff3_feats['product'] = ''
            if not 'note' in gff3_feats:
                gff3_feats['note'] = ''
            print >> gff3_fh, "%s\tgenbankfile\texon\t%d\t%d\t.\t%s\t0\tID=%s.exon.%s%s;Parent=%s.transcript.%s%s" % (formatogff(rec.id), gff3_feats['start'], gff3_feats['end'], gff3_feats['strand'], formatogff(gff3id), input_number, featid, formatogff(gff3id), input_number, featid)
            print >> gff3_fh, "%s\tgenbankfile\tmRNA\t%d\t%d\t.\t%s\t.\tID=%s.transcript.%s%s;Parent=%s.gene.%s%s%s%s;Note=%s;Name=%s.transcript.%s%s;%s" % (formatogff(rec.id), gff3_feats['start'], gff3_feats['end'], gff3_feats['strand'], formatogff(gff3id), input_number, featid, formatogff(gff3id), input_number, featid, gff3_feats['EC'], gff3_featsCDS, formatogff(gff3_feats['product']),  formatogff(gff3id), input_number, featid, gff3_feats['note'])
            print >> gff3_fh, "%s\tgenbankfile\tCDS\t%d\t%d\t.\t%s\t0\tID=%s.CDS.%s%s;Parent=%s.transcript.%s%s" % (formatogff(rec.id), gff3_feats['start'], gff3_feats['end'], gff3_feats['strand'], formatogff(gff3id), input_number, featid, formatogff(gff3id), input_number, featid)
            gff3_featsCDS = ''

        if  f.type == 'tRNA':
            if not 'product' in gff3_feats:
                gff3_feats['product'] = ''
                if 'note' in gff3_feats:
                    gff3_feats['product'] = gff3_feats['note']

            print >> gff3_fh, "%s\tgenbankfile\texon\t%d\t%d\t.\t%s\t0\tID=%s.exon.%s%s" % (formatogff(rec.id), gff3_feats['start'], gff3_feats['end'], gff3_feats['strand'], formatogff(gff3id), input_number, featid)

            print >> gff3_fh, "%s\tgenbankfile\ttRNA\t%d\t%d\t.\t%s\t.\tID=%s.transcript.%s%s;Note=%s;description=%s" % (formatogff(rec.id), gff3_feats['start'], gff3_feats['end'], gff3_feats['strand'], formatogff(gff3id), input_number, featid, formatogff(gff3_feats['product']), formatogff(gff3_feats['product']))
            print >> gff3_fh, "%s\tgenbankfile\tCDS\t%d\t%d\t.\t%s\t0\tID=%s.CDS.%s%s" % (formatogff(rec.id), gff3_feats['start'], gff3_feats['end'], gff3_feats['strand'],formatogff(gff3id), input_number, featid)
        if  f.type == 'rRNA':
            print >> gff3_fh, "%s\tgenbankfile\texon\t%d\t%d\t.\t%s\t0\tID=%s.exon.%s%s" % (formatogff(rec.id), gff3_feats['start'], gff3_feats['end'], gff3_feats['strand'], formatogff(gff3id), input_number, featid)
            print >> gff3_fh, "%s\tgenbankfile\trRNA\t%d\t%d\t.\t%s\t.\tID=%s.transcript.%s%s;Note=%s;description=%s" % (formatogff(rec.id), gff3_feats['start'], gff3_feats['end'], gff3_feats['strand'], formatogff(gff3id), input_number, featid, formatogff(gff3_feats['product']), formatogff(gff3_feats['product']))
            print >> gff3_fh, "%s\tgenbankfile\tCDS\t%d\t%d\t.\t%s\t0\tID=%s.CDS.%s%s" % (formatogff(rec.id), gff3_feats['start'], gff3_feats['end'], gff3_feats['strand'],formatogff(gff3id), input_number, featid)



feature_fh.close()
fasta_fh.close()
gff3_fh.close()
