from pymongo import MongoClient
import pandas as pd
from Bio.PDB.PDBParser import PDBParser
from Bio.SeqUtils import seq1
from glob import glob
from SNDG import init_log
from SNDG.BioMongo.Model.SeqCollection import SeqCollection
from SNDG.BioMongo.Model.Protein import Protein
from SNDG.BioMongo.Model.Alignment import AlnLine
from SNDG.BioMongo.Model.ResidueAln import ResidueAln
from SNDG.BioMongo.Model.Structure import ModeledStructure, Chain, Molecule, StructureQuality, \
    ResidueSet
from SNDG.BioMongo.Process.StructureAnotator import StructureAnotator
from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
from SNDG.BioMongo.Process.Index import StructuromeIndexer
from SNDG.WebServices.Offtargeting import  Offtargeting
from SNDG.BioMongo.Process.Importer import load_pathways,import_prop_blast,build_statistics,index_seq_collection

import pymongo
from bson import ObjectId
from tqdm import tqdm
import Bio.SearchIO as bpsio
import logging
import json

init_log()
_log = logging.getLogger(__name__)
import os

# with open('/data/databases/pdb/manually_curated/compound_type.csv') as handle:
#     compound_type = {x.replace('"', "").strip().split(",")[0]: x.replace('"', "").strip().split(",")[1] if
#     x.replace('"', "").strip().split(",")[1] else "?" for x in handle.readlines()}
#
#
# def get_compound_type(residue):
#     return compound_type[residue.get_resname().strip()] if residue.get_resname().strip() in compound_type else "?"
#
#
# parser = PDBParser(PERMISSIVE=1, QUIET=1)
# basepath = "/data/projects/structsN315/results/"
# # model_files = glob(basepath + "*/*/*/*.pdb")
#
# organism = "SaureusN315"
# seq_col_id = ObjectId("5b2800b1be737e35a6dd9b8a")
#
# port = 27018
# mdb = BioMongoDB("tdr",port=port)
# db = pymongo.MongoClient(port=port).pdb
# name = "SaureusN315"


# with tqdm(model_files) as pbar:
#     for model_file in pbar:
#         pbar.set_description("processing %s" % model_file)
#         model_name = model_file.split("/")[-2]
#         seq_name = model_file.split("/")[-3]
#         aln = [hit[0] for hit in list(bpsio.read(basepath + "/" + seq_name + "/profile_search.xml", "blast-xml")) if
#                hit.id == model_name.split(seq_name + "_")[1]][0]
#
#         with open(model_file + ".json") as h:
#             assessments = json.load(h)
#         pockets = []
#
#         prot = list(Protein.objects(organism=organism, gene=seq_name))
#         if len(prot) == 0:
#             _log.warn("Not found: " + seq_name)
#         try:
#             structure = parser.get_structure(model_name, model_file)
#             with open(model_file) as h:
#                 pdblines = [l.strip() for l in h.readlines() if l.startswith("REMARK")]
#             strdoc = ModeledStructure(name=model_name, seq_collection_id=seq_col_id, pipeline="sndg",
#                                       organism=organism)
#
#             squality = StructureQuality(name="ga341", value=float(
#                 [l.split("GA341 score: ")[1] for l in pdblines if "GA341 score: " in l][0]))
#             strdoc.qualities.append(squality)
#             squality = StructureQuality(name="dope", value=float(
#                 [l.split("DOPE score: ")[1] for l in pdblines if "DOPE score: " in l][0]))
#             strdoc.qualities.append(squality)
#             squality = StructureQuality(name="qmean", value=assessments["QMEAN4_norm"])
#             strdoc.qualities.append(squality)
#             squality = StructureQuality(name="zqmean", value=assessments["QMEAN4_zscore"])
#             strdoc.qualities.append(squality)
#
#             model = structure[0]
#             for chain in model:
#                 chaindoc = Chain(name=chain.id)
#                 strdoc.chains.append(chaindoc)
#
#                 aln_h = str(aln.aln[1].seq)
#                 aln_q = str(aln.aln[0].seq)
#                 prot_start = aln.query_start
#                 prot_end = aln.query_end
#                 pdb_txt = [l for l in pdblines if "TEMPLATE: " in l][0].split()
#                 pdb_start = int(pdb_txt[4].split(":")[0])
#                 pdb_end = int(pdb_txt[6].split(":")[0])
#
#                 hit_start = aln.hit_start
#                 hit_end = aln.hit_end
#
#                 template_name = aln.aln[1].id
#
#                 aln_query = AlnLine(name=seq_name, seq=aln_q, start=prot_start, end=prot_end)
#                 aln_hit = AlnLine(name=template_name, seq=aln_h, start=hit_start, end=hit_end)
#                 template = ResidueAln(aln_query=aln_query, aln_hit=aln_hit,
#                                       query_res_start=1,
#                                       query_res_end=len(chain),
#                                       hit_res_start=pdb_start,
#                                       hit_res_end=pdb_end)
#                 strdoc.templates.append(template)
#
#                 for residue in chain:
#                     molecule = Molecule(resid=residue.id[1],
#                                         chain=chain.id,
#                                         compound=seq1(residue.get_resname()),
#                                         compound_type=get_compound_type(residue))
#                     if molecule.compound_type == 'RESIDUE':
#                         chaindoc.residues.append(molecule)
#                     else:
#                         if not [x for x in strdoc.ligands if
#                                 (
#                                         x.compound_type == molecule.compound_type) and molecule.compound_type == 'SOLVENT']:
#                             strdoc.ligands.append(molecule)
#
#             pockets_json = model_file + ".pockets.json"
#             if os.path.exists(pockets_json):
#                 rss = StructureAnotator.pocket_residue_set(pockets_json, model.get_atoms())
#                 strdoc.pockets = rss
#             strdoc.save()
#         except Exception as ex:
#             _log.error(ex)





# with tqdm(glob("/data/organismos/SaureusN315/estructura/sndg/pockets/*.json")) as pbar:
#     for pockets_json in pbar:
#         model = pockets_json.split("/")[-1].split(".")[0]
#         pbar.set_description("processing %s" % model)
#         strdoc = ModeledStructure.objects(organism="SaureusN315",name=model).get()
#         structure = parser.get_structure(model, "/data/organismos/SaureusN315/estructura/sndg/modelos/" + model + ".pdb")
#         with open(pockets_json) as h:
#             pocket_data = json.load(h)
#         if os.path.exists(pockets_json)  and pocket_data:
#             rss = StructureAnotator.pocket_residue_set(pockets_json, structure[0].get_atoms())
#             strdoc.pockets = rss
#             strdoc.save()
#
#

# sa = StructureAnotator(basepath, struct_path=lambda wd, modeldoc: glob( "/".join(
#     [wd, modeldoc.templates[0].aln_query.name, modeldoc.templates[0].aln_query.name, modeldoc.name, "*.pdb"]))[0])
# total = sa.total(db, organism, {})
#
# with tqdm(sa.iterator(db, organism, {}), total=total) as pbar:
#     for model in pbar:
#         pbar.set_description(model.name)
#
#         template = model.templates[0]
#         protein = Protein.objects(organism=organism, alias=template.aln_query.name).get()
#
#         sa.annotate_model(model, protein.domains())
#         model.save()
#
#
# si = StructuromeIndexer(SeqCollection.objects(name=name).get())
# si.build_index()

# load_pathways("SaureusN315", "/data/organismos/SaureusN315/annotation/pwtools/smallmolecules.sbml", mdb.db,
#               "/data/organismos/SaureusN315/annotation/pwtools/",
#               gregexp="\(([\-\w\.]+)\)", filter_file="allfilters_con_c.dat", gene_map="/data/organismos/SaureusN315/annotation/pwtools/gene_map.pkl")



# Offtargeting.offtargets("/data/organismos/SaureusN315/annotation/proteins.fasta",
#                         "/data/organismos/" + name + "/annotation/offtarget/",
#                         cpus=2  )
# import_prop_blast(mdb.db, name, "hit_in_deg",
#                   "/data/organismos/" + name + "/annotation/offtarget/degaa-p.tbl",
#                   "table", "Hit in DEG database",
#                   value_fn=lambda x: x.identity > 70,
#                   default_value="True",
#                   no_hit_value=False, choices=["True", "False"], type="value", defaultOperation="equal")

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

# from SNDG.BioMongo.Process.BioCyc2Mongo import BioCyc
# biocyc = BioCyc(mdb.db)
# biocyc.user = BioMongoDB.demo
# biocyc.pre_build_index(SeqCollection.objects(name=name).get())
#
# build_statistics(mdb.db,name)