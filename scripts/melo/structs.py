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
import pymongo
from bson import ObjectId
from tqdm import tqdm
import Bio.SearchIO as bpsio
import logging
import json

init_log()
_log = logging.getLogger(__name__)
import os

with open('/data/databases/pdb/manually_curated/compound_type.csv') as handle:
    compound_type = {x.replace('"', "").strip().split(",")[0]: x.replace('"', "").strip().split(",")[1] if
    x.replace('"', "").strip().split(",")[1] else "?" for x in handle.readlines()}


def get_compound_type(residue):
    return compound_type[residue.get_resname().strip()] if residue.get_resname().strip() in compound_type else "?"


parser = PDBParser(PERMISSIVE=1, QUIET=1)
organism = "Linf"
basepath = "/data/projects/infantum/estructuras/"
model_files = [] # glob(basepath + "*/*/*/*.pdb")

# model_files = [
#     basepath + "Minc3s00001g00003/Minc3s00001g00003/Minc3s00001g00003_3uaf_A_21_137/Minc3s00001g00003.B99990001.pdb"]
models_count = len(model_files)

seq_col_id = ObjectId("5b2800b1be737e35a6dd9b8a")

BioMongoDB("tdr")
db = pymongo.MongoClient().pdb

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
#             pockets_json = model_file + "pockets.json"
#             if os.path.exists(pockets_json):
#                 rss = StructureAnotator.pocket_residue_set(pockets_json, model.get_atoms())
#                 strdoc.pockets = rss
#             strdoc.save()
#         except Exception as ex:
#             _log.error(ex)

sa = StructureAnotator(basepath, struct_path=lambda wd, modeldoc: glob( "/".join(
    [wd, modeldoc.templates[0].aln_query.name, modeldoc.templates[0].aln_query.name, modeldoc.name, "*.pdb"]))[0])
total = sa.total(db, organism, {})



with tqdm(sa.iterator(db, organism, {}), total=total) as pbar:
    for model in pbar:
        pbar.set_description(model.name)


        template = model.templates[0]
        try:
            protein = Protein.objects(organism=organism, alias=template.aln_query.name).get()
        except :
            print template.aln_query.name + " does not exists"
            continue
        sa.annotate_model(model, protein.domains())
        model.save()


