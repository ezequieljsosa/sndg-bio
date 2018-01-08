#!/usr/bin/python
import os
import sys
#unip_goa_path = sys.argv[1] if len(sys.argv) > 1 else  "gp_association.goa_uniprot"

unip_goa_path = "goa_uniprot_all.gpa"


goa = open(unip_goa_path)
goaheader = open("gp_association.goa_uniprot.000", "w")
ciclo=1
goaline = goa.next()
print "	generando gp_association.goa_uniprot.000"
while ciclo:
	if goaline.split("\t")[0] != "UniProtKB": 
		goaheader.write(goaline)
		goaline = goa.next()
	else:
		prefix =  "".join(list(goaline.split("\t")[1])[0:3])
		print  "	generando gp_association.goa_uniprot.%s" % prefix
		goapart =  open("gp_association.goa_uniprot.%s" % prefix, "w")
		while "".join(list(goaline.split("\t")[1])[0:3]) == prefix:
			goapart.write(goaline)
			try: 
				goaline = goa.next()
			except: 
				ciclo=0
				break
		goapart.close()
goa.close()
goaheader.close()

for i in range(65,91):
	for j in range(0,10):
		for k in range(65,91):
			if not os.path.exists("gp_association.goa_uniprot.%s%s%s" % (chr(i), str(j), chr(k)) ):
				new_gp = open("gp_association.goa_uniprot.%s%s%s" % (chr(i), str(j), chr(k) ),"w")
				print "Archivo mock gp_association.goa_uniprot.%s%s%s" % (chr(i), str(j), chr(k) )
				new_gp.close()

sys.exit(0)
