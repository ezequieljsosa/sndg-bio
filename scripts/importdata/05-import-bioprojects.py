import logging
import traceback

from Bio import Entrez
from tqdm import tqdm

from SNDG import init_log, mkdir
from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
from SNDG.BioMongo.Process.Importer import from_ref_seq, update_proteins, create_proteome
from SNDG.BioMongo.Process.Taxon import tax_db
from SNDG.WebServices.NCBI import ExternalAssembly, mysql_db, Submitter, BioProject
from peewee import MySQLDatabase
from SNDG.Sequence.ProteinAnnotator import ProteinAnnotator, Mapping
from SNDG.BioMongo.Process.Index import index_seq_collection, build_statistics

Entrez.email = "ezejajaja@hotmail.com"
_log = logging.getLogger(__name__)


mysql_db.initialize(MySQLDatabase('sndg', user='root', passwd="mito"))

# submitters = Submitter.select().where((Submitter.source == "ncbi") & ( Submitter.rejected == False ))

BioProjectSubmitter = BioProject.submitters.get_through_model()

BioProjectSubmitter.drop_table()
BioProject.drop_table()

BioProject.create_table()
BioProjectSubmitter.create_table()




esearch = Entrez.read(Entrez.esearch(db="bioproject", term="argentina", retmax=10000))["IdList"]
def to_utf8(string):
    try:
        return string.decode('utf-8')
    except UnicodeError:
        return string.encode("utf-8").decode("utf-8")

for pid in tqdm(esearch):
    bioproject = Entrez.read(Entrez.esummary(db="bioproject", id=pid))["DocumentSummarySet"]["DocumentSummary"][0]
    submitters = []
    for x in bioproject["Submitter_Organization_List"]:
        x = to_utf8(x)
        if list(Submitter.select().where(Submitter.name == x)):
            s = Submitter.select().where(Submitter.name == x).get()
            submitters.append(s)
        else:
            s = Submitter(name=x, source="ncbi")
            s.save()
            submitters.append(s)



    bp = BioProject(
        identifier=to_utf8(bioproject["Project_Acc"]),
        name=to_utf8(bioproject["Project_Title"]),
        description=to_utf8(bioproject["Project_Description"]),
        material=to_utf8(bioproject["Project_Target_Capture"]),
        scope=to_utf8(bioproject["Project_Target_Material"])

    )
    try:
        bp.save()

    except Exception as ex:
        bp.description = ""
        bp.save()
    finally:
        bp.submitters.add(submitters)
        bp.save()


