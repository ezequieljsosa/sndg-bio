import pymongo

from Bio import Entrez
from tqdm import tqdm
from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
from SNDG.WebServices.NCBI import ExternalAssembly, mysql_db, Submitter, BioProject, NCBI, AssemblySubmitters
from peewee import MySQLDatabase

mysql_db.initialize(MySQLDatabase('sndg', user='root', passwd="mito"))

mdb = BioMongoDB("saureus")
# local_submitter
# submitters
assemblies = list(ExternalAssembly.select())
for assembly in tqdm(assemblies):
    if not mdb.db.bioprojects.count({"ncbi_accession": assembly.assembly_accession}):
        mdb.db.sequence_collection.update()
        mdb.db.proteins.update()
        mdb.db.contig_collection.update()

