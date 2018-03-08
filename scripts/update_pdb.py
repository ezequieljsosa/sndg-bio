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
import json
import logging
import os
import sys
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

import pandas as pd
from Bia import init_log
from Bia.Model.collections import SeqCollection
from BiaStructure.IO import PDBsIterator
from BiaStructure.Model.structure import ExperimentalStructure, Chain, Molecule, \
    ResidueSet
from BiaStructure.Programs.fpocket import fpocket_properties_map
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import CaPPBuilder
from pymongo.mongo_client import MongoClient
from tqdm import tqdm

from SNDG.Structure.CompoundTypes import compound_type, get_compound_type
from SNDG import mkdir

fpocket_properties_map_inv = {v: k for k, v in fpocket_properties_map.items()}

init_log()
_log = logging.getLogger(__name__)

entries_columns = ["IDCODE", "HEADER", "ACCESSIONDATE", "COMPOUND", "SOURCE", "AUTHORS", "RESOLUTION", "EXPERIMENT"]
pdb_entries_df = pd.read_table('/data/databases/pdb/entries.idx', skiprows=[0, 1, 2], sep='\t', names=entries_columns)


def complete_pdb_attrs(pdb, struct_doc):
    # TODO QMEAN SASA FOLDX DSSP BINDINGDB DRUGDATABANK!!
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


def complete_pockets(pdb, strdoc):
    pdb_file = "/data/databases/pdb/divided/%s/pdb%s.ent" % (pdb[1:3], pdb)
    pockets_json = "/data/databases/pdb/pockets/%s/%s.ent.json" % (pdb[1:3], pdb)
    mkdir("/data/databases/pdb/pockets/%s/" % (pdb[1:3]))

    if not os.path.exists(pockets_json) or os.path.getsize(pockets_json) < 10:
        r = FPocket(pdb_file).hunt_pockets()
        r.save(pockets_json)
        r.result.delete_dir()

    rs = StructureAnotator.pocket_residue_set(pockets_json)
    if rs.druggability >= 0.2:
        strdoc.pockets.append(rs)


def procesar_pdb(pdb, parser, ppb, collection):
    pdb_file = "/data/databases/pdb/divided/%s/pdb%s.ent" % (pdb[1:3], pdb)
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
            strdoc = ExperimentalStructure(name=pdb, seq_collection_id=collection, seq_collection_name="pdb")

            model = structure[0]
            for chain in model:
                chaindoc = Chain(name=chain.id,
                                 segments=[[y.id[1] for y in list(x)] for x in ppb.build_peptides(chain)])
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

            complete_pockets(pdb, strdoc)

            complete_pdb_attrs(pdb, strdoc)
            strdoc.save()
    except Exception as ex:
        with open("pdb_load_errors.txt", "a") as handle:
            handle.write(pdb + "|" + str(ex) + "\n")


if __name__ == "__main__":
    argv = sys.argv

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]

    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose", dest="verbose", action="count",
                        help="set verbosity level [default: %(default)s]")

    parser.add_argument("-db", "--database_name", default='tdr')
    parser.add_argument("-host", "--db_host", default='127.0.0.1')
    parser.add_argument('-V', '--version', action='version', version=program_version_message)

    args = parser.parse_args()
    verbose = args.verbose

    if verbose > 0:
        print("Verbose mode on")

    db = MongoClient(args.db_host)[args.database_name]
    col_name = "pdb"

    from mongoengine import connect

    connect(args.database_name, host=args.db_host)
    collection = SeqCollection.objects(name=col_name)
    if len(collection):
        collection = collection.get()
    else:
        collection = SeqCollection(name=col_name, description="Protein Data Bank", organism="?")
        collection.save()

    parser = PDBParser(PERMISSIVE=1, QUIET=1)
    ppb = CaPPBuilder()

    procesados_sin_pocket = [x["name"] for x in
                             db.structures.find({"seq_collection_name": "pdb", "pockets.0": {"$exists": 0}},
                                                {"name": 1})]
    procesados = [x["name"] for x in db.structures.find({"seq_collection_name": "pdb"}, {"name": 1})]
    pdbs = list(PDBsIterator())
    for (pdb, pdb_file) in tqdm(pdbs):

        if pdb in procesados:
            if pdb in procesados_sin_pocket:
                q = ExperimentalStructure.objects(seq_collection_name="pdb", name=pdb).no_cache()
                if q:
                    strdoc = q.get()
                    try:
                        complete_pockets(pdb, strdoc)
                    except Exception as ex:
                        _log.warm("ERROR " + pdb + " : " + str(ex))
                        continue

                    if strdoc.pockets:
                        strdoc.save()
        else:
            procesar_pdb(pdb, parser, ppb, collection)

    print "OK!"
