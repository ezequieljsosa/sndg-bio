import pymongo

from Bio import Entrez
from tqdm import tqdm
from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
from SNDG.WebServices.NCBI import ExternalAssembly, mysql_db, Submitter, BioProject, NCBI, AssemblySubmitters
from peewee import MySQLDatabase

mysql_db.initialize(MySQLDatabase('sndg', user='root', passwd="mito"))
from SNDG.BioMongo.Process.KeywordIndexer import  KeywordIndexer
mdb = BioMongoDB("saureus")
# local_submitter
# submitters
# assemblies = list(ExternalAssembly.select())
# for assembly in tqdm(assemblies):
#     is_local = len([x for x in assembly.submitters if not x.submitter.rejected])
#     mdb.db.sequence_collection.update({"ncbi_assembly": assembly.assembly_accession},
#                                       {"$set":{"submitters":[x.submitter.name for x in assembly.submitters]}})
#     if is_local:
#         add_local_kw = {"$addToSet": {"keywords": "local_submitter"}}
#         mdb.db.sequence_collection.update({"ncbi_assembly": assembly.assembly_accession},
#                                           add_local_kw)
#
#         record = mdb.db.sequence_collection.find_one({"ncbi_assembly": assembly.assembly_accession}, {"name": 1})
#         if record:
#             mdb.db.proteins.update({"organism": record["name"]}, add_local_kw)
#             mdb.db.contig_collection.update({"organism": record["name"]}, add_local_kw)
ki = KeywordIndexer()
bioprojects = list(BioProject.select())
for biop in tqdm(bioprojects):
    is_local = len([x for x in biop.submitters if not x.rejected])
    keywords = ki.extract_keywords(biop.description)
    if is_local:
        keywords.append("local_submitter")

    print mdb.db.bioprojects.update({"accession": biop.accession},
                                        {"$set":{"keywords":keywords}})