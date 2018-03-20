# or x in SeqCollection.objects():
#    ...:     if x.tax:
#    ...:             db = "/data/databases/deg/deg-e-13.3/degaa-e.dat"
#    ...:             if x.tax.superkingdom == "procaryotae":
#    ...:                 db = "/data/databases/deg/deg-p-13.3/degaa-p.dat"
#    ...:             adir = "/data/organismos/" + x.name + "/analysis/deg/"
#    ...:             if not os.path.exists(adir):
#    ...:                 os.makedirs(adir)
#    ...:             h = open(adir + "cmd.txt","w")
#    ...:
#    ...:             cmd = ("blastp -query /data/organismos/" + x.name + "/anotacion/proteins.fasta -db %s -evalue 0.00001 -outfmt 6 -out " + adir  + "deg_hits.tbl -max_hsps 1 -qcov_hsp_perc 0.8 -max_target_seqs
#    ...: 1") % db
#    ...:             print cmd
#    ...:             h.write(cmd)
#    ...:             h.close()
#    ...:             sp.call(cmd, shell=True)


# In [9]: for x in SeqCollection.objects():
#    ...:     if x.tax:
#    ...:         print x.name
#    ...:         adir = "/data/organismos/" + x.name + "/anotacion/uniprot/"
#    ...:         prots = "/data/organismos/" + x.name + "/anotacion/proteins.fasta"
#    ...:         os.chdir(adir)
#    ...:         sp.call("makeblastdb -in sp.fasta -dbtype prot",shell=True)
#    ...:         sp.call("makeblastdb -in tr.fasta -dbtype prot",shell=True)
#    ...:         templ = "blastp -db {base}.fasta -query {prots} -evalue 0.00001 -max_hsps 1 -qcov_hsp_perc 0.8 -outfmt 5 > {base}_blast.xml"
#    ...:         sp.call( templ.format(base="sp",prots=prots)   ,shell=True)
#    ...:         sp.call( templ.format(base="tr",prots=prots)   ,shell=True)