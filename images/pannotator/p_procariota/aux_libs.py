#!/usr/bin/python
"""
Created on Fri Jun 19 19:05:31 2015

@author: gburguener
"""

import re
import os
import sys
import time
import datetime
import threading
import getopt
from random import randint
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.NCBIStandalone import BlastParser
from Bio.Blast.NCBIStandalone import Iterator
from Bio.Blast import NCBIXML
from commands import getoutput
from math import floor
import multiprocessing

def cpu_count(): 
	return multiprocessing.cpu_count()

def tiempo_transcurrido(begin):
	now = datetime.datetime.now()
	segundosnow = now.day*24*60*60 + now.hour*60*60 + now.minute*60 + now.second
	segundosbegin = begin.day*24*60*60 + begin.hour*60*60 + begin.minute*60 + begin.second
	dias = floor((segundosnow - segundosbegin) / (24*60*60))
	horas = floor((segundosnow - segundosbegin - dias*24*60*60)/(60*60))
	minutos = floor((segundosnow - segundosbegin - dias*24*60*60 - horas*60*60)/60)
	segundos = segundosnow - segundosbegin - dias*24*60*60 - horas*60*60 - minutos*60
	return [dias, horas, minutos, segundos]


class BackgroundBlast(object):
	def __init__(self, query, database, output="blastout", threads="4", evalue="1e-5", program="blastp"):
		self.query = query
		self.database = database
		self.output = output
		self.threads = str(threads)
		self.evalue = str(evalue)
		self.program = program
		
		self.background = threading.Thread(target=self.blast, args=())
		self.background.daemon = True

	def blast(self):
		getoutput("blastall -p %s -i '%s' -d '%s' -e %s -a %s -v 20 -b 20 -o '%s' " % (self.program, self.query, self.database, self.evalue, self.threads, self.output))
		return    0
		
		
	def run(self):
		self.background.start()
