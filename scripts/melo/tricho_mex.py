import logging
from SNDG import Struct, init_log, mkdir

init_log()
from glob import glob
import pymongo
from SNDG.WebServices.NCBI import NCBI
from SNDG.WebServices.Offtarget import Offtarget
import Bio.SeqIO as bpio
from tqdm import tqdm
import json




from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB

from SNDG.BioMongo.Process.Importer import from_ref_seq, update_proteins, import_prop_blast,common_annotations
from SNDG.BioMongo.Process.BioDocFactory import BioDocFactory
from SNDG.BioMongo.Model.Protein import Protein
from SNDG.Network.KEGG import Kegg
from SNDG.BioMongo.Process.Importer import _common_annotations, _protein_iter, import_kegg_annotation, \
    index_seq_collection, build_statistics, load_pathways,common_annotations
from BCBio import GFF

from SNDG.BioMongo.Model.Structure import ModeledStructure, Molecule, ResidueAln, SimpleAlignment, StructureQuality, \
    ExperimentalStructure,Chain
from SNDG.BioMongo.Model.Alignment import AlnLine
import os
from SNDG.BioMongo.Process.StructureAnotator import StructureAnotator
import Bio.SearchIO as bpsio

mdb = BioMongoDB("tdr", port=27017)


name = "tatro"
organism = name
org = "Trichoderma atroviride"
ann_path = "/data/organismos/tatro/annotation/corrected.gb"
# mdb.delete_seq_collection(name)
# from_ref_seq(name, ann_path,  cpus=6)
# common_annotations(name, "/data/organismos/tatro/annotation/", cpu=6, remove_tmp=False)
# mdb.protein_fasta("/data/organismos/" + name + "/annotation/proteins.faa", name)
from SNDG.Annotation.EMapper import EMapper
# em = EMapper()
# em.read_file("proteins.")
#update_proteins("/tmp/" + name + "/", "/data/organismos/" + name + "/annotation/proteins.faa", name, 1003200, db_init=mysqldb)
#
#
# Offtargeting.offtargets("/data/organismos/" + name + "/annotation/proteins.faa",
#                         "/data/organismos/" + name + "/annotation/offtarget/",
#                         offtarget_dbs=[  "/data/databases/deg/degaa-p.dat",
#                                          "/data/databases/human/gencode.v17.pc_translations.fa",
#                                          "/data/databases/human/gut_microbiota.fasta"]
#                         )
# import_prop_blast(mdb.db, name, "hit_in_deg",
#                   "/data/organismos/" + name + "/annotation/offtarget/degaa-p.tbl",
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


# index_seq_collection(mdb.db,name,pathways=False,go=False,keywords=False,ec=False,organism_idx=False,structure=True)
# load_pathways(name, "/data/organismos/Ainsu/annotation/pwtools/small_molecule.sbml", mdb.db,
#               "/data/organismos/Ainsu/annotation/pathways",
#               gregexp="\(([\-\w\.]+)\)", filter_file="allfilters_con_c.dat")
# index_seq_collection(mdb.db,name,pathways=False,go=False,keywords=True,ec=False,organism_idx=True,structure=False)
# build_statistics(mdb.db,name)
dst_db = pymongo.MongoClient(port=27018).tdr
mdb.copy_genome(name,name,dst_db)