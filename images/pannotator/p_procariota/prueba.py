#!/usr/bin/python

import re
from random import randint
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from commands import getoutput



map = open("/data/uniprot/genesymbolmap")
genome_iterator=SeqIO.parse(map ,"fasta")

new = open("/data/uniprot/genesymbol", "w")

while 1:
	try:
		line = map.next()
	except:
		break
	id = line.split("|")[1]
	GN =  line.split("GN=")[-1].split()[0]
	PE =  line.split("PE=")[-1].split()[0]
	
	print >> new, id, GN, PE 