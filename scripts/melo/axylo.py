import logging
from SNDG import Struct, init_log, mkdir

init_log()
import glob
import pymongo
from SNDG.WebServices.NCBI import NCBI
from SNDG.WebServices.Offtargeting import Offtargeting
import Bio.SeqIO as bpio
from tqdm import tqdm

logging.getLogger("peewee").setLevel(logging.WARN)
from peewee import MySQLDatabase
from SNDG.BioMongo.Process.Taxon import tax_db
from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
from SNDG.Sequence.ProteinAnnotator import ProteinAnnotator, Mapping
from SNDG.BioMongo.Process.Importer import from_ref_seq, update_proteins, import_prop_blast
from SNDG.BioMongo.Process.BioDocFactory import BioDocFactory
from SNDG.BioMongo.Model.Protein import Protein
from SNDG.Network.KEGG import Kegg
from SNDG.BioMongo.Process.Importer import _common_annotations, _protein_iter, import_kegg_annotation, \
    index_seq_collection, build_statistics, load_pathways
from BCBio import GFF
from SNDG.BioMongo.Process.Taxon import Tax

tax_db.initialize(MySQLDatabase('bioseqdb', user='root', passwd="mito"))
mdb = BioMongoDB("tdr",port=27018)
mysqldb = ProteinAnnotator.connect_to_db(database="unipmap", user="root", password="mito")

name = "Axylo"
org = "Achromobacter xylosoxidans"
ann_path = "/data/organismos/Axylo/annotation/GCF_001457475.1_NCTC10807_genomic.gbff"

# from_ref_seq(name, ann_path,  tax=85698,      cpus=2)

# mdb.protein_fasta("/data/organismos/Axylo/annotation/proteins.faa", name)
# update_proteins("/tmp/" + name + "/", "/data/organismos/Axylo/annotation/proteins.faa", name, 85698, db_init=mysqldb)


# Offtargeting.offtargets("/data/organismos/Axylo/annotation/proteins.faa",
#                         "/data/organismos/" + name + "/annotation/offtarget/",
#                         offtarget_dbs=[  "/data/databases/deg/degaa-e.dat",
#                                          "/data/databases/human/gencode.v17.pc_translations.fa",
#                                          "/data/databases/human/gut_microbiota.fasta"]
#                         )
# import_prop_blast(mdb.db, name, "hit_in_deg",
#                   "/data/organismos/" + name + "/annotation/offtarget/degaa-e.tbl",
#                   "table", "Hit in DEG database",
#                   value_fn=lambda x: x.identity > 70,
#                   default_value=True,
#                   no_hit_value=False, choices=["True", "False"], type="value", defaultOperation="equal")
#
# import_prop_blast(mdb.db, name, "human_offtarget",
#                   "/data/organismos/" + name + "/annotation/offtarget/gencode.tbl",
#                   "table", "Human offtarget score (1 - best hit identity)",
#                   value_fn=lambda x: 1 - (x.identity * 1.0 / 100),
#                   default_value=0.4,
#                   no_hit_value=1)
# import_prop_blast(mdb.db, name, "gut_microbiota_offtarget",
#                   "/data/organismos/" + name + "/annotation/offtarget/gut_microbiota.tbl",
#                   "table", "Gut microbiota offtarget score (1 - best hit identity)",
#                   value_fn=lambda x: 1 - (x.identity * 1.0 / 100),
#                   default_value=0.4,
#                   no_hit_value=1)

# kegg = Kegg(ko_list_path="/data/databases/kegg/ko.txt",
#             brittle_pw_path="/data/databases/kegg/ko/",
#             kgmls_dir="/data/databases/kegg/ko/")
# kegg.init()
# ilex_data = "/data/organismos/" + name + "/query.ko"
# kegg.read_annotation(ilex_data)
#
# import_kegg_annotation(mdb.db,name,kegg)



# index_seq_collection(mdb.db,name,pathways=False,go=False,keywords=False,ec=False,organism_idx=False,structure=True)
# load_pathways(name, "/data/organismos/Axylo/annotation/pwtools/all.sbml", mdb.db, "/data/organismos/Axylo/annotation/pathways",
#               gregexp="\(([\-\w\.]+)\)", filter_file="allfilters_con_c.dat")
# index_seq_collection(mdb.db,name,pathways=True,go=True,keywords=True,ec=True,organism_idx=True,structure=False)
# build_statistics(mdb.db,name)
