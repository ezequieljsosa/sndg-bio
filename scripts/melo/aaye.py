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
    index_seq_collection, build_statistics, load_pathways
from BCBio import GFF

from SNDG.BioMongo.Model.Structure import ModeledStructure, Molecule, ResidueAln, SimpleAlignment, StructureQuality, \
    ExperimentalStructure,Chain
from SNDG.BioMongo.Model.Alignment import AlnLine
import os
from SNDG.BioMongo.Process.StructureAnotator import StructureAnotator
import Bio.SearchIO as bpsio

mdb = BioMongoDB("tdr", port=27017)


name = "Aaye"
organism = name
org = "Acinetobacter baumannii AYE"
ann_path = "/mnt/data/workspace/aaye/GCF_000069245.1_ASM6924v1_genomic.gbff"

# from_ref_seq(name, ann_path,  cpus=3)
# common_annotations(name, "./", cpu=4, remove_tmp=False)
# mdb.protein_fasta("/data/organismos/" + name + "/annotation/proteins.faa", name)
# from SNDG.Annotation.EMapper import EMapper
# em = EMapper()
# em.read_file("proteins.")
# update_proteins("/tmp/" + name + "/", "/data/organismos/" + name + "/annotation/proteins.faa", name, )
mdb = BioMongoDB("tdr")
# mdb.load_from_emapper("Absce","./abse.emapper.annotation")
#mdb.load_from_interpro(name,"/mnt/data/data/organismos/Aaye/annotation/interproscan-5.46-81.0/proteins.faa.gff3")


# go2ec = {}
# for line in open("/mnt/data/data/organismos/Aaye/annotation/ec2go.txt"):
#     if line.startswith("EC:"):
#         ec = line.split()[0].lower()
#         go = line.strip().split()[-1].lower()
#         ec = ec + "".join([".-" for _ in range(3-ec.count("."))])
#         go2ec[go] = ec
# for p in mdb.db.proteins.find({"organism":name,"ontologies":{'$regex':"^go"}}):
#     for go in p["ontologies"]:
#         if "go:" in go and go in go2ec:
#             mdb.db.proteins.update({"_id":p["_id"]},{"$addToSet":{"ontologies":go2ec[go]}})


# Offtarget.offtargets(f"/mnt/data/data/organismos/{name}/annotation/proteins.faa",
#                      f"/mnt/data/data/organismos/{name}/annotation/offtarget/deg.tbl",
#                      "/data/databases/deg/degaa-p.dat", cpus=3)
# Offtarget.offtargets(f"/mnt/data/data/organismos/{name}/annotation/proteins.faa",
#                      f"/mnt/data/data/organismos/{name}/annotation/offtarget/human.tbl",
#                      "/data/databases/human/human_uniprot100.fa.gz", cpus=3)

# import_prop_blast(mdb.db, name, "hit_in_deg",
#                   "/data/organismos/" + name + "/annotation/offtarget/deg.tbl",
#                   "table", "Hit in DEG database",
#                   value_fn=lambda x: x.identity > 70,
#                   default_value=True,
#                   no_hit_value=False, choices=["True", "False"], type="value", defaultOperation="equal")
#
# import_prop_blast(mdb.db, name, "human_offtarget",
#                   "/data/organismos/" + name + "/annotation/offtarget/human.tbl",
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

#add column names
# mdb.load_metadata( name,"/mnt/data/data/organismos/Aaye/annotation/gut_microbiome.tbl")


# from pymongo import MongoClient
# db = MongoClient().pdb
# sa = StructureAnotator(basepath, struct_path=lambda wd, modeldoc: glob("/".join(
#     [wd, modeldoc.templates[0].aln_query.name, modeldoc.templates[0].aln_query.name, modeldoc.name, "*.pdb"]))[0])
# total = sa.total(db, organism, {})
#
# with tqdm(sa.iterator(db, organism, {}), total=total) as pbar:
#     for model in pbar:
#         pbar.set_description(model.name)
#
#         template = model. templates[0]
#         protein = Protein.objects(organism=organism, alias=template.aln_query.name).get()
#         sa.annotate_model(model, protein.domains())
#         model.save()


# index_seq_collection(mdb.db,name,pathways=False,go=False,keywords=False,ec=False,organism_idx=False,structure=True)
# load_pathways(name, "/data/organismos/Aaye/annotation/pathways/small.sbml", mdb.db,
#               "/data/organismos/Aaye/annotation/pathways_data/",
#               gregexp="\(([\-\w\.]+)\)", filter_file="allfilters_con_c.dat")
# index_seq_collection(mdb.db,name,pathways=True,go=True,keywords=True,ec=True,organism_idx=True,structure=False)
# build_statistics(mdb.db,name)


mdb.copy_genome("Aaye",dst_db=pymongo.MongoClient(port=27018).tdr)