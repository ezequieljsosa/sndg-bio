import pymongo

from Bio import Entrez
from tqdm import tqdm
from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
from SNDG.WebServices.NCBI import ExternalAssembly, mysql_db, Submitter, BioProject, NCBI, AssemblySubmitters
from peewee import MySQLDatabase

mysql_db.initialize(MySQLDatabase('sndg', user='root', passwd="mito"))

mdb = BioMongoDB("saureus")

bioprojects = list(BioProject.select())
for biop in tqdm(bioprojects):
    if not mdb.db.bioprojects.count({"accession": biop.accession}):
        biop_doc = {
            "name": biop.name,
            "identifier": biop.identifier,
            "accession": biop.accession,
            "material": biop.material,
            "scope": biop.scope,
            "description": biop.description,
            "created": biop.created,
            "submitters": [s.name for s in biop.submitters],
            "assemblies": [{"accession": a.assembly_accession, "name": a.name} for a in biop.assemblies],
        }
        mdb.db.bioprojects.insert(biop_doc)
