#!/usr/bin/env python

map = open("./uniprotheaders")
new = open("./genesymbolmap", "w")

while 1:
	try:
		line = map.next()
	except:
		break

	id = line.split("|")[1]
	GN =  line.split("GN=")[-1].split()[0]
	PE =  line.split("PE=")[-1].split()[0]
	
	print id
	print >> new, id, GN, PE 
