import os
import logging
import traceback

from SNDG.WebServices import download_file
from tqdm import tqdm
from SNDG import execute,mkdir
from SNDG.WebServices.NCBI import NCBI

from Bio import Entrez

_log = logging.getLogger(__name__)




class Offtargeting(object):

    offtarget_p = [  "/data/databases/deg/degaa-p.dat",
                     "/data/databases/human/gencode.v17.pc_translations.fa",
                     "/data/databases/human/gut_microbiota.fasta"]

    @staticmethod
    def download_deg(dst="/data/databases/deg/"):
        for x in ["p", "e", "a"]:
            filename = "deg-" + x + "-15.2"

            download_file("http://tubic.tju.edu.cn/deg/download/" + filename + ".zip",
                          dst +  filename + ".zip", ovewrite=True)
            execute("unzip -o  " + dst  + filename + ".zip" + " -d " + dst)
            os.remove(dst + filename + ".zip")
            execute("makeblastdb -dbtype prot -in " + dst + "degaa-" + x + ".dat")

    @staticmethod
    def download_human_prots(dst="/data/databases/human/"):
        download_file("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_17/gencode.v17.pc_translations.fa.gz",
                      dst +  "gencode.v17.pc_translations.fa.gz", ovewrite=True)
        execute("gunzip --directory > " + dst + "gencode.v17.pc_translations.fa")

    @staticmethod
    def create_human_microbiome(dst="/data/databases/human/",
                                accession_list="/data/databases/human/accessions_gut_microbiota_2017.txt"):
        Entrez.email = "anemail@hotmail.com"
        dst_accs = dst + "gut_microbiota_assemblies/"
        mkdir(dst_accs)
        out_file = dst + "gut_microbiota_raw.fasta"
        errors = []
        if os.path.exists(out_file):
            os.remove(out_file)
        with open(accession_list) as h:
            for accession in h.readlines():
                try:
                    accession = accession.strip().split(".")[0]
                    proteome_path = NCBI.download_assembly(accession, dst_accs,dtype="protein.faa.gz")
                    execute("cat " + proteome_path + " >> " + out_file)
                    os.remove(proteome_path)
                except:
                    errors.append(accession)
        final_file = dst + "gut_microbiota.fasta"
        execute("cd-hit -c 0.95 -o " + final_file + " -i " + out_file)
        execute("makeblastdb -dbtype prot -in " + final_file)
        print errors
        return out_file



    # @staticmethod
    # def offtargets(proteome,dst_resutls,offtarget_dbs=):
    #     for db in []

if __name__ == "__main__":
    from SNDG import init_log
    from SNDG.WebServices import PROXIES
    PROXIES["ftp_proxy"] = "http://proxy.fcen.uba.ar:8080"
    init_log()
    #Offtargeting.download_deg()
    #Offtargeting.download_human_prots()
    Offtargeting.create_human_microbiome()

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
