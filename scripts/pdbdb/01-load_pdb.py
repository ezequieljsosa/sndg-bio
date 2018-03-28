'''
Created on Nov 2, 2017

@author: eze
'''

import logging
from argparse import ArgumentParser, RawDescriptionHelpFormatter

import math
from Bio.PDB.PDBParser import PDBParser

from tqdm import tqdm

from SNDG import init_log
from SNDG.Structure.PDBdb import *
from SNDG.Structure.PDBs import PDBs

init_log(rootloglevel=logging.INFO)
_log = logging.getLogger(__name__)

if __name__ == "__main__":

    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)

    parser.add_argument("-p", "--dbpass", required=True)
    parser.add_argument("-i", "--pdb_dir", default="/data/databases/pdb/")
    parser.add_argument("-db", "--dbname", default="pdbdb")
    parser.add_argument("-u", "--dbuser", default="root")

    args = parser.parse_args()
    from peewee import MySQLDatabase

    mysql_db = MySQLDatabase(args.dbname, user=args.dbuser, password=args.dbpass)
    mysql_db.close()
    sqldb.initialize(mysql_db)

    pdb_utils = PDBs(pdb_dir=args.pdb_dir)
    df = pdb_utils.entries_df()
    pdbs = list(pdb_utils)
    with tqdm(pdbs) as pbar:
        for (code, pdb_path) in pbar:
            mysql_db.connect(reuse_if_open=True)
            pbar.set_description(code)
            try:
                entry = df[df.IDCODE == code.upper()].iloc[0]
            except IndexError:
                continue

            pdb_model = PDB(code=code, experiment=str(entry.EXPERIMENT))

            try:
                resolution = float(entry.RESOLUTION)
            except:
                if not math.isnan(resolution):
                    pdb_model.resolution = resolution
            pdb_model.save()

            p = PDBParser(PERMISSIVE=True, QUIET=True)
            try:
                for chain in p.get_structure(code, pdb_path).get_chains():
                    with sqldb.atomic():
                        residues = []
                        for residue in chain.get_residues():
                            if residue.resname != "HOH":
                                residue_model = {"pdb": pdb_model, "chain": chain.id, "resid": residue.id[1],
                                                 "icode": residue.id[2],
                                                 "type": "RES" if not residue.id[0].strip() else residue.id[0],
                                                 "resname": residue.resname, "disordered": residue.is_disordered()}
                                residues.append(residue_model)
                        Residue.insert_many(residues).execute()
                    with sqldb.atomic():
                        residues = list(Residue.select().where((Residue.pdb == pdb_model) &
                                                               (Residue.chain == chain.id)))
                        atoms = []
                        for idx, residue in enumerate(chain.get_residues()):
                            if residue.resname != "HOH":
                                for atom in residue.get_atoms():
                                    residue_model = residues[idx]
                                    assert atom.parent.id[1] == residue_model.resid

                                    atm = {"residue":residue_model, "serial":atom.serial_number, "name":atom.id,
                                         "x":float(atom.coord[0]), "y":float(atom.coord[1]), "z":float(atom.coord[2]),
                                         "occupancy":float(atom.occupancy), "bfactor":float(atom.bfactor),
                                         "element":atom.element}
                                    atoms.append(atm)
                        Atom.insert_many(atoms).execute()



            except Exception as ex:
                _log.error(code + ": " + str(ex))
            finally:
                mysql_db.close()
