import os
import logging
import traceback

from SNDG.WebServices import download_file
from tqdm import tqdm
from SNDG import execute, mkdir
from SNDG.WebServices.NCBI import NCBI
import multiprocessing

from Bio import Entrez

_log = logging.getLogger(__name__)

from collections import defaultdict
class Offtargeting(object):
    offtarget_p = ["/data/databases/deg/degaa-p.dat",
                   "/data/databases/human/gencode.v17.pc_translations.fa",
                   "/data/databases/human/gut_microbiota.fasta"]

    @staticmethod
    def download_deg(dst="/data/databases/deg/"):
        for x in ["p", "e", "a"]:
            filename = "deg-" + x + "-15.2"

            download_file("http://tubic.tju.edu.cn/deg/download/" + filename + ".zip",
                          dst + filename + ".zip", ovewrite=True)
            execute("unzip -o  " + dst + filename + ".zip" + " -d " + dst)
            os.remove(dst + filename + ".zip")
            execute("makeblastdb -dbtype prot -in " + dst + "degaa-" + x + ".dat")

    @staticmethod
    def download_human_prots(dst="/data/databases/human/"):
        download_file(
            "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_17/gencode.v17.pc_translations.fa.gz",
            dst + "gencode.v17.pc_translations.fa.gz", ovewrite=True)
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
                    proteome_path = NCBI.download_assembly(accession, dst_accs, dtype="protein.faa.gz")
                    execute("cat " + proteome_path + " >> " + out_file)
                    os.remove(proteome_path)
                except:
                    errors.append(accession)
        final_file = dst + "gut_microbiota.fasta"
        execute("cd-hit -c 0.95 -o " + final_file + " -i " + out_file)
        execute("makeblastdb -dbtype prot -in " + final_file)
        print errors
        return out_file



    @staticmethod
    def count_organism_from_microbiome_blast(tbl_blast_result_path, microbiome_fasta,identity_threshold=0.4, out_tbl=None):
        prot_org_map = {}
        with open(microbiome_fasta) as h:
            for line in h:
                if line.startswith(">"):
                    seqid = line.split()[0].strip().replace(">","")
                    org = line.replace("[[", "[").split("[")[1].strip()[:-1]
                    prot_org_map[seqid] = org

        query_orgs = defaultdict(lambda: [])
        with open(tbl_blast_result_path) as h:
            for l in list(h)[1:]:
                query, hit,identity = l.split()[:3]
                identity = float(identity) / 100.0
                if identity_threshold <= identity:
                    query_orgs[query].append(prot_org_map[hit])
        for query, hits in query_orgs.items():
            query_orgs[query] = set(hits)
        if out_tbl:
            with open(out_tbl,"w") as h:
                for query, hits in query_orgs.items():
                    h.write("\t".join([query, str(len(hits)), ";".join(hits)]) + "\n")
        return query_orgs

    @staticmethod
    def offtargets(proteome, dst_resutls, offtarget_dbs=offtarget_p, cpus=multiprocessing.cpu_count()):
        mkdir(dst_resutls)
        pname = proteome.split("/")[-1].split(".")[0]
        results = []
        for db in offtarget_dbs:
            dbname = db.split("/")[-1].split(".")[0]
            out = dst_resutls + dbname + ".tbl"
            execute(
                "blastp -evalue 1e-5 -max_hsps 1 -outfmt 6 -max_target_seqs 1 -db {db} -query {query} -out {out} -num_threads {cpus}",
                db=db, query=proteome, out=out, cpus=cpus)
            results.append(out)
        return results


if __name__ == "__main__":
    from SNDG import init_log
    from SNDG.WebServices import PROXIES

    PROXIES["ftp_proxy"] = "http://proxy.fcen.uba.ar:8080"
    init_log()
    # Offtargeting.download_deg()
    # Offtargeting.download_human_prots()
    # Offtargeting.create_human_microbiome()

    Offtargeting.offtargets("/data/organismos/GCF_001624625.1/annotation/genoma.fasta",
                            "/data/organismos/GCF_001624625.1/annotation/offtarget/"
                            )
