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
mdb = BioMongoDB("tdr")
mysqldb = ProteinAnnotator.connect_to_db(database="unipmap", user="root", password="mito")

name = "Linf"
org = "Leishmania infantum LinJ36"
ann_path = "/data/projects/infantum/GCA_900500625.1_LINF_genomic.gbff"

# from_ref_seq(name, ann_path, seqs=None, tax=5671, tmp_dir=None,
#              extract_annotation_feature=lambda feature: feature.sub_features[
#                  0] if feature.type == "gene" and hasattr(feature,"sub_features") and len(feature.sub_features) else feature,
#              accept_protein_feature=lambda f: ((f.type == "CDS") and ("translation" in f.qualifiers)),
#              extract_sequence=lambda c, f: f.extract(c).seq.translate(),
#              cpus=4)

# mdb.protein_fasta("/tmp/" + name + "/genome.fasta", name)
# update_proteins("/tmp/" + name + "/", "/tmp/" + name + "/genome.fasta", name, 5658, db_init=mysqldb)


# Offtargeting.offtargets("/tmp/" + name + "/genome.fasta",
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
# load_pathways(name, "/data/organismos/Linf/annotation/pwtools/allreactions.sbml", mdb.db, "/data/organismos/Linf/annotation/pwtools/",
#               gregexp="\(([\-\w\.]+)\)", filter_file="allfilters_con_c.dat")
# index_seq_collection(mdb.db,name,pathways=False,go=True,keywords=True,ec=True,organism_idx=True,structure=True)
# build_statistics(mdb.db,name)
