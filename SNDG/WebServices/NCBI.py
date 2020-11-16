#!/usr/bin/env python

'''
Created on Feb 21, 2017

@author: eze
'''

import datetime
import logging
import os
import subprocess as sp
import sys
import requests
from tqdm import tqdm

from Bio import Entrez

Entrez.email = 'A.N.Other@example.com'
from SNDG import execute, init_log
from SNDG.WebServices import download_file

ftp_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all"

_log = logging.getLogger(__name__)


class NCBI(object):
    f_mRNA = "mRNA"
    f_CDS = "CDS"
    ftypes = ["rRNA", "ncRNA", f_mRNA, "gene", "exon", f_CDS, "rRNA", "tRNA", "tmRNA"]

    # https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly
    dbs = ["bioproject", "biosample", "dbvar", "sra", "assembly"]  # "nuccore","protein",
    # "genome"
    dbs_con_submitter = ["assembly"]  # "bioproject", "biosample",

    # "pubmed", ?
    # "gene", ?

    def __init__(self):
        '''
        Constructor
        '''

        self.resource_handler = {
            "assembly": NCBIAssembly(), "bioproject": NCBIProject(), "biosample": NCBIBiosample(),
            "gene": NCBIGene(), "genome": NCBIGenome(), "dbvar": None, "sra": NCBIReads(), "pubmed": NCBIPubmed(),
            "protein": NCBIProtein(), "nuccore": NCBINucleotide()

        }

    @staticmethod
    def download_assembly(assembly_accession, dst_dir, dtype="genomic.gbff.gz", force=False):
        # assembly_name, last_assembly_accession = NCBI.assembly_name_from_acc(assembly_accession)
        assembly_accession_no_ver = assembly_accession if assembly_accession[-2] != "." else assembly_accession[:-2]

        # https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/158/435/
        url = "/".join([ftp_url
                           , assembly_accession_no_ver[0:3], assembly_accession_no_ver[4:7]
                           , assembly_accession_no_ver[7:10], assembly_accession_no_ver[10:13]
                        # , last_assembly_accession + "_" + assembly_name.replace(" ", "_").replace("#", "_")
                        # , last_assembly_accession + "_" + assembly_name.replace(" ", "_").replace("#",
                        #                                                                           "_") + "_" + dtype
                        ]) + "/"
        r = requests.get(url)
        download_url = ""
        acc =""
        if r.status_code == 200:
            accessions = [x.split("</")[0].replace("/","") for x in r.text.split(">") if x.startswith(assembly_accession_no_ver)]
            # GCF_000158435.(1)_ASM15843v1/
            accessions = sorted(accessions, key=lambda x: int(x.split("_")[1][-1]))
            if accessions:
                acc = accessions[-1]
                download_url = f'{url}{acc}/{acc}_{dtype}'

        if not download_url:
            err = f"{assembly_accession} not found at {url}"
            _log.error(err)
            raise FileNotFoundError(err)

        assert acc

        out_file = f'{dst_dir}/{"_".join(acc.split("_")[:2]) }.{dtype}'
        if force or not os.path.exists(out_file):
            download_file(download_url, out_file, ovewrite=force)
        else:
            _log.debug(f'{out_file} exists')
        # execute("gunzip -c  " + out_file + " > " +  out_file[:-3])
        return out_file

    def download(self, accesion, db, dst, dstformat):
        cmd = 'wget -O %s "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=%s&report=%s&sort=&id=%s&from=begin&to=end"'
        params = (dst, db, dstformat, accesion)

        execute(cmd % params, shell=True)

    @staticmethod
    def assembly_name_from_acc(assembly_accession):
        # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term=GCA_000009645.1
        esearch = Entrez.read(Entrez.esearch(db="assembly", term=assembly_accession))

        if esearch["IdList"]:
            assembly_id = esearch["IdList"][0]
            #             https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id=31148
            summary = Entrez.esummary(db="assembly", id=assembly_id)
            summary = Entrez.read(summary, validate=False)
            ata = summary["DocumentSummarySet"]["DocumentSummary"][0]
            last_ver = max([int(data[x].split(".")[-1]) for x in
                            ['LatestAccession', 'LastMajorReleaseAccession', 'AssemblyAccession'] if
                            x in data and data[x]])
            assembly_accession = assembly_accession if "." not in assembly_accession else ".".join(
                assembly_accession.split(".")[:-1])
            return (data["AssemblyName"], assembly_accession + "." + str(last_ver))
        else:
            raise FileNotFoundError(assembly_accession)


def org_name(org_complete_name):
    arr = org_complete_name.replace("sequence", "").replace("genome", "").replace(",", "").replace(".", "").replace(" ",
                                                                                                                    "").replace(
        ":", "").split(" ")
    if len(arr) == 1:
        name = org_complete_name
    elif len(arr) == 2:
        name = arr[0][0:4] + arr[1][0:4]
    else:
        name = arr[0][0:4] + arr[1][0:4]
        if "subsp" in arr:
            if "strain" in arr:
                name += "".join(arr[arr.index("subsp") + 1: arr.index("strain")][:2])
            else:
                name += "".join(arr[arr.index("subsp") + 1: arr.index("subsp") + 3])
        if "strain" in arr:
            name += "".join(arr[arr.index("strain") + 1: arr.index("strain") + 3])
        else:
            name += arr[-1] if arr[-1] not in arr[0:2] else ""
    return name


if __name__ == "__main__":

    import argparse
    import Bio.SeqIO as bpio
    import os
    import fileinput
    import subprocess as sp
    import traceback

    parser = argparse.ArgumentParser(description='Download NCBI assemblies')
    parser.add_argument('-a', '--accession', required=True)
    parser.add_argument('-o', '--output', help="output_directory", default="./")
    parser.add_argument('--force', action="store_true")

    args = parser.parse_args()

    if args.accession == "-":
        for accession in fileinput.input("-"):
            accession = accession.strip()
            sys.stderr.write(f'{accession}\n')
            try:
                downloaded_file = NCBI.download_assembly(accession, dst_dir=args.output)
                if not os.path.exists(downloaded_file):
                    sys.stderr.write(f"error downloading {accession}\n")
                elif os.path.getsize(downloaded_file) < 10:
                    os.remove(downloaded_file)
                    sys.stderr.write(f"error downloading {accession}\n")
                else:
                    sys.stdout.write(downloaded_file)
                    sys.stdout.write("\n")
            except sp.CalledProcessError:
                traceback.print_exc(file=sys.stderr)
    else:
        downloaded_file = NCBI.download_assembly(args.accession, dst_dir=args.output, force=args.force)
        if not os.path.exists(downloaded_file) or (os.path.getsize(downloaded_file) < 10):
            sys.stderr.write(f"error downloading {args.accession}")
        else:
            sys.stdout.write(downloaded_file)
            sys.stdout.write("\n")

    # logger = logging.getLogger('peewee')
    # logger.setLevel(logging.INFO)
    # init_log()
    # # Submitter.create_table()
    # # ExternalResource.create_table()
    # ExternalAssembly.create_table()
    # AssemblySubmitters.create_table()
    #
    # Entrez.email = "ezejajaja@hotmail.com"
    # # NCBI().actualizar_submitters(NCBI().obtener_submitters())
    #
    # submitters = Submitter.select().where((Submitter.source == "ncbi") & (Submitter.rejected == False))
    # NCBI().update_assemblies(submitters)
    #
    # sys.exit()
    # nombres = []
    # # s = Submitter.select(Submitter.name == "sndg").get()
    #
    # # usehistory="y"
    # buscar = True
    # init = 0
    # while buscar:
    #     retstart = init
    #     pepe = Entrez.read(Entrez.esearch(db="assembly", retmax="200",
    #                                       term='"representatives"[Filter]AND"latest refseq"[Filter]AND"bacteria"[Filter]AND"complete genome"[Assembly Level]'
    #                                       , usehistory="Y", retstart=str(retstart)))
    #     if pepe["IdList"]:
    #         juan = Entrez.read(Entrez.esummary(db="assembly", id=",".join(pepe["IdList"])), validate=False)
    #         for assembly in juan["DocumentSummarySet"]["DocumentSummary"]:
    #             # genome = str(Entrez.read(Entrez.elink(dbfrom="assembly", id=assembly.attributes["uid"], linkname="assembly_genome")) [0]['IdList'][0])
    #             genome = str(assembly["SpeciesName"])
    #             assert assembly['LastMajorReleaseAccession'].startswith("GCF")
    #             if not ExternalAssembly.select().where(
    #                     ExternalAssembly.assembly_accession == str(assembly["AssemblyAccession"])).count():
    #                 ExternalAssembly(assembly_accession=str(assembly["AssemblyAccession"]),
    #                                  assembly_name=str(assembly["AssemblyName"]),
    #                                  genome=genome,
    #                                  submitter_id=74,
    #                                  name=genome + " ref assembly",
    #                                  type="assembly",
    #                                  identifier=str(assembly.attributes["uid"])
    #                                  ).save()
    #                 # "representatives"[Filter]  "latest refseq"[Filter] "bacteria"[Filter] "complete genome"[Assembly Level] "chromosome"[Assembly Level]
    #         init = init + 200
    #     else:
    #         buscar = False

#         workdir = "/data/projects/ncbi_dump/" 
#         os.chdir(workdir)
#         for x in ExternalAssembly.select():
#             aworkdir= str(x.id)  
#             file_name =  x.gbk_file_name().replace(".gz","")                        
#             if not os.path.exists(aworkdir):                
#                 os.makedirs( aworkdir )
#                 x.download_gbk(aworkdir)                
#             result = aworkdir + "/" + file_name
#             if os.path.exists(result):
#                 org =  org_name(list(bpio.parse(result,"gb"))[0].description)
#                 print org
#                 nombres.append(org)
#             else:
#                 print "no exists: " + x.identifier
#         print len(nombres)
#         print len(set(nombres))

# from SNDG.BioMongo.Process.Taxon import Tax, tax_db
# tax_db.initialize(MySQLDatabase('bioseqdb', user='root', passwd="mito"))
# tax = Tax.getTax(1872703)
# for col in db.sequence_collection.find({"tax":{"$exists":0}, "ncbi_assembly":{"$exists":1} },{"name":1,"ncbi_assembly":1,"_id":0}):
#     esearch = Entrez.read(Entrez.esearch(db="assembly", term= col["ncbi_assembly"] ))["IdList"][0]
#     summary = Entrez.esummary(db="assembly", id=esearch)
#     summary = Entrez.read(summary, validate=False)
#     juan = dict(summary["DocumentSummarySet"]["DocumentSummary"][0])
#     #pepe["AssemblyStatus"]
#     tax = Tax.getTax(str(juan["Taxid"]))
#     tax = {
#         "tid" : float(tax.ncbi_taxon_id),
#                          "superkingdom" : [y for y in [x for x in Tax.parents(tax) if x.node_rank == "superkingdom"][0].names
#                                            if y.name_class == "scientific name"  ][0].name  ,
#         "name" : [x for x in tax.names if x.name_class == "scientific name"][0].name
#         }
#     print tax
#     print db.sequence_collection.update({"name":col["name"]},{"$set":{"tax":tax}})
