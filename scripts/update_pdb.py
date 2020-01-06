#!/usr/bin/python
# encoding: utf-8
'''
scripts.load_bia_proteome -- shortdesc

scripts.load_bia_proteome is a description

It defines classes_and_methods

@author:     Ezequiel Sosa
@copyright:  2015 BIA. All rights reserved.
@license:    license
@contact:    user_email

'''
import logging
import os
import sys
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

from Bio.PDB.PDBParser import PDBParser
from pymongo.mongo_client import MongoClient
from tqdm import tqdm

from SNDG import init_log
from SNDG import mkdir
from SNDG.BioMongo.Model.Structure import ExperimentalStructure, Chain, Molecule
from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
from SNDG.BioMongo.Process.StructureAnotator import StructureAnotator
from SNDG.Structure.CompoundTypes import get_compound_type
from SNDG.Structure.FPocket import FPocket
from SNDG.Structure.PDBs import PDBs

init_log(rootloglevel=logging.INFO)
_log = logging.getLogger(__name__)

pdbUtils = PDBs()


def complete_pdb_attrs(pdb, struct_doc):
    # TODO QMEAN SASA FOLDX DSSP BINDINGDB DRUGDATABANK!!
    pdb_entries_df = pdbUtils.entries_df()
    row = pdb_entries_df[pdb_entries_df.IDCODE == pdb.upper()]
    if len(row) > 0:
        row = row.iloc[0]
        struct_doc.description = row["COMPOUND"]
        try:
            if row["SOURCE"].strip() != "":
                struct_doc.organism = row["SOURCE"].strip()
        except:
            pass
        if row["RESOLUTION"].strip():
            try:
                struct_doc.resolution = float(row["RESOLUTION"].strip())
            except:
                pass
        if row["EXPERIMENT"].strip() != "":
            struct_doc.experiment = row["EXPERIMENT"].strip()


def complete_pockets(pdb, strdoc, structure):
    pdb_file = pdbUtils.pdb_path(pdb)
    pockets_json = pdbUtils.pdb_pockets_path(pdb)
    mkdir("/data/databases/pdb/pockets/%s/" % (pdb[1:3]))

    if not os.path.exists(pockets_json) or os.path.getsize(pockets_json) < 10:
        r = FPocket(pdb_file).hunt_pockets()
        r.save(pockets_json)
        r.delete_dir()

    if os.path.exists(pockets_json):
        strdoc.pockets = StructureAnotator.pocket_residue_set(pockets_json, structure.get_atoms())


def procesar_pdb(pdb, parser):
    pdb_file = pdbUtils.pdb_path(pdb)
    if not os.path.exists(pdb_file):
        with open("/tmp/pdb_load_errors.txt", "a") as handle:
            handle.write(pdb + "|NO EXISTE: " + pdb_file + " \n")
    try:
        structure = parser.get_structure(pdb, pdb_file)
        models = list(structure)
        if len(models) == 0:
            with open("/tmp/pdb_load_errors.txt", "a") as handle:
                handle.write("No tiene modelos: " + pdb_file + " \n")
        else:
            strdoc = ExperimentalStructure(name=pdb, seq_collection_name="pdb")

            model = structure[0]
            for chain in model:
                chaindoc = Chain(name=chain.id,
                                 segments=[[y.id[1] for y in list(x)] for x in [chain]])
                strdoc.chains.append(chaindoc)

                for residue in chain:
                    res_id = residue.id[1]
                    molecule = Molecule(resid=res_id,
                                        chain=chain.id,
                                        compound=residue.get_resname(),
                                        compound_type=get_compound_type(residue))
                    if get_compound_type(residue) == 'RESIDUE':
                        chaindoc.residues.append(molecule)

                    else:
                        molecule.compound_type = get_compound_type(residue)
                        #                         if not [x for x in strdoc.ligands if (x.compound_type == molecule.compound_type) and  (x.compound_type == 'SOLVENT')]:
                        if molecule.compound_type != 'SOLVENT':
                            strdoc.ligands.append(molecule)
            try:
                complete_pockets(pdb, strdoc, structure)
            except:
                pass

            complete_pdb_attrs(pdb, strdoc)
            strdoc.save()
    except Exception as ex:
        with open("pdb_load_errors.txt", "a") as handle:
            handle.write(pdb + "|" + str(ex) + "\n")


if __name__ == "__main__":
    argv = sys.argv

    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose", dest="verbose", action="count",
                        help="set verbosity level [default: %(default)s]")

    parser.add_argument("-host", "--db_host", default='127.0.0.1')
    BioMongoDB("sndg")
    args = parser.parse_args()

    db = MongoClient(args.db_host)["pdb"]
    col_name = "pdb"

    """
    collection = SeqCollection.objects(name=col_name)
    if len(collection):
        collection = collection.get()
    else:
        collection = SeqCollection(name=col_name, description="Protein Data Bank", organism="?")
        collection.save()
    """
    parser = PDBParser(PERMISSIVE=1, QUIET=1)

    procesados_sin_pocket = {x["name"]: 1 for x in
                             db.structures.find({"seq_collection_name": "pdb", "pockets.0": {"$exists": 0}},
                                                {"name": 1})}
    procesados = {x["name"]: 1 for x in db.structures.find({"seq_collection_name": "pdb"}, {"name": 1})}
    pdbs = list(PDBs())
    for (pdb, pdb_file) in tqdm(pdbs):

        if pdb in procesados:
            if pdb in procesados_sin_pocket:
                q = ExperimentalStructure.objects(seq_collection_name="pdb", name=pdb).no_cache()
                if q:
                    strdoc = q.get()
                    try:
                        structure = parser.get_structure(pdb, pdbUtils.pdb_path(pdb))
                        procesar_pdb(pdb, strdoc, structure)
                    except Exception as ex:
                        _log.warn("ERROR " + pdb + " : " + str(ex))
                        continue
                    if strdoc.pockets:
                        strdoc.save()
        else:
            procesar_pdb(pdb, parser)

    print ("OK!")
