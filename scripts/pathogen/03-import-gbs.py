from tqdm import tqdm

from SNDG import mkdir
from SNDG.WebServices.NCBI import NCBI
from Bio import Entrez, SeqIO
import Bio.SeqIO as bpio
import multiprocessing
import subprocess as sp

if __name__ == '__main__':
    Entrez.email = "ezequieljsosa@gmail.com"
    query = '"pathogen"[Properties] AND ("Metazoa"[Organism] OR "Viridiplantae"[Organism] OR "Fungi"[Organism] OR "Eukaryota"[Organism] NOT "Metazoa"[Organism] NOT "Fungi"[Organism] NOT "Viridiplantae"[Organism] OR "Bacteria"[Organism] OR txid1224[Orgn] OR "Archaea"[Organism])'
    genomesList = Entrez.read(Entrez.esearch(db="genome", term=query, idtype="acc", retmax=10000))

    genomes = Entrez.read(Entrez.esummary(db="genome", id=",".join(genomesList["IdList"])), validate=False)
    with tqdm(genomes) as pbar1:
        for genome in pbar1:
            pbar1.set_description(genome["Organism_Name"])
            query = '"%s"[Organism] AND ((latest[filter] OR "latest refseq"[filter]) AND all[filter] NOT anomalous[filter])'
            query += " AND " + '("complete genome"[filter] OR "chromosome level"[filter] OR "scaffold level"[filter])'
            query = query % genome["Organism_Name"]

            try:
                assembliesList = Entrez.read(Entrez.esearch(db="assembly", term=query, idtype="acc", retmax=100000))
                for assembly_id in assembliesList["IdList"]:
                    resource = Entrez.read(Entrez.esummary(db="assembly", id=assembly_id), validate=False)
                    data = resource["DocumentSummarySet"]["DocumentSummary"][0]
                    acc = data["AssemblyAccession"]
                    asName = data["AssemblyName"]
                    basedir = "/tmp/" + genome["Organism_Name"] + "/" + acc + "/"

                    asspath = "/".join([acc[0:3], acc[4:7], acc[7:10], acc[10:13],
                                        acc + "_" + asName.replace(" ", "_").replace("#", "_")])

                    cmd = 'rsync --recursive --include="*_genomic.gbff.gz" --exclude="*"   rsync://ftp.ncbi.nlm.nih.gov/genomes/all/' + asspath + '/ "' + basedir + '"'
                    sp.call(cmd, shell=True)
                    # def save_acc():
                    #     sequencesList = Entrez.read(Entrez.esearch(
                    #         db="nuccore", term=str(acc) + "[Assembly]", idtype="acc", retmax=100000))["IdList"]
                    #     basedir = "/tmp/" + genome["Organism_Name"] + "/"
                    #     seq_ids = [str(x) for x in sequencesList]
                    #     mkdir(basedir)
                    #     def seq_iterator():
                    #
                    #             handle = Entrez.efetch(db="nuccore", id=",".join(seq_ids), rettype="gbwithparts", retmode="text")
                    #             for seq in SeqIO.parse(handle, "genbank"):
                    #                 yield seq
                    #     pepe= list(seq_iterator())
                    #     bpio.write(pepe, basedir + acc + ".gb", "gb")
                    #
                    #
                    # p = multiprocessing.Process(target=save_acc)
                    # p.start()
                    # p.join(180)
                    # if p.is_alive():
                    #     p.terminate()
                    #     p.join()
                    #     raise multiprocessing.TimeoutError()


            except:
                with open("errors.log", "a") as h:
                    h.write(str(genome["Organism_Name"]) + "\n")
