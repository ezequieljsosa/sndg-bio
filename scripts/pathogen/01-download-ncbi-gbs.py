from tqdm import tqdm
import os
from SNDG import mkdir
from SNDG.WebServices.NCBI import NCBI
from Bio import Entrez, SeqIO
import Bio.SeqIO as bpio
import multiprocessing
import subprocess as sp
from glob import glob
import urllib2
import httplib

FNULL = open(os.devnull, 'w')

if __name__ == '__main__':
    Entrez.email = "ezequieljsosa@gmail.com"
    query = '"pathogen"[Properties] AND ("Metazoa"[Organism] OR "Viridiplantae"[Organism] OR "Fungi"[Organism] OR "Eukaryota"[Organism] NOT "Metazoa"[Organism] NOT "Fungi"[Organism] NOT "Viridiplantae"[Organism] OR "Bacteria"[Organism] OR txid1224[Orgn] OR "Archaea"[Organism])'
    genomesList = Entrez.read(Entrez.esearch(db="genome", term=query, idtype="acc", retmax=10000))
    processed = {}
    if os.path.exists("processed.txt"):
        processed = {x.strip(): 1 for x in open("processed.txt").readlines()}
    genomes = Entrez.read(Entrez.esummary(db="genome", id=",".join(genomesList["IdList"])), validate=False)
    with tqdm(genomes) as pbar1:
        for genome in pbar1:
            pbar1.set_description(genome["Organism_Name"])
            query = '"%s"[Organism] AND ((latest[filter] OR "latest refseq"[filter]) AND all[filter] NOT anomalous[filter])'
            query += " AND " + '("complete genome"[filter] OR "chromosome level"[filter] OR "scaffold level"[filter])'
            query = query % genome["Organism_Name"]
            try:
                assembliesList = Entrez.read(Entrez.esearch(db="assembly", term=query, idtype="acc", retmax=100000))
            except  (urllib2.URLError, httplib.BadStatusLine):
                continue
            for assembly_id in tqdm(assembliesList["IdList"]):
                if assembly_id in processed:
                    continue
                try:
                    resource = Entrez.read(Entrez.esummary(db="assembly", id=assembly_id), validate=False)
                    data = resource["DocumentSummarySet"]["DocumentSummary"][0]
                    acc = data["AssemblyAccession"]
                    asName = data["AssemblyName"]
                    basedir = "/data/bio/patho2/" + genome["Organism_Name"] + "/"

                    asspath = "/".join([acc[0:3], acc[4:7], acc[7:10], acc[10:13],
                                        acc + "_" + asName.replace(" ", "_").replace("#", "_")])
                    if glob(basedir + acc + "*"):
                        with open("processed.txt", "a") as h:
                            h.write(str(assembly_id) + "\n")
                        continue

                    cmd = 'rsync --recursive --include="*_genomic.gbff.gz" --exclude="*"   rsync://ftp.ncbi.nlm.nih.gov/genomes/all/' + asspath + '/ "' + basedir + '"'
                    sp.call(cmd, shell=True, stdout=FNULL)
                    with open("processed.txt", "a") as h:
                        h.write(str(assembly_id) + "\n")

                except  (urllib2.URLError, httplib.BadStatusLine):
                    with open("errors.log", "a") as h:
                        h.write(str(genome["Organism_Name"]) + " " + str(assembly_id) + "\n")
