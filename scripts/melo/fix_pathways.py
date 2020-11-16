import logging
from SNDG import Struct, init_log, mkdir

init_log()
import glob
import pymongo

from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
from SNDG.BioMongo.Process.Importer import _common_annotations, _protein_iter, import_kegg_annotation, \
    index_seq_collection, build_statistics, load_pathways

mdb = BioMongoDB("tdr")

load_pathways("MeloI", "/data/organismos/MeloI/annotation/pwtools/allreactions.sbml", mdb.db,
              "/data/organismos/MeloI/annotation/pwtools/",
              gregexp="\(([\-\w\.]+)\)", filter_file="allfilters_con_c.dat")

load_pathways("SaureusN315", "/data/organismos/SaureusN315/annotation/pwtools/smallmolecules.sbml", mdb.db,
              "/data/organismos/SaureusN315/annotation/pwtools/",
              gene_map="/data/organismos/SaureusN315/annotation/pwtools/gene_map.pkl",
              gregexp="\(([\-\w\.]+)\)", filter_file="allfilters_con_c.dat")

load_pathways("Linf", "/data/organismos/Linf/annotation/pwtools/allreactions.sbml", mdb.db,
              "/data/organismos/Linf/annotation/pwtools/",
              gregexp="\(([\-\w\.]+)\)", filter_file="allfilters_con_c.dat")