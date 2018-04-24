'''
Created on Nov 2, 2017

@author: eze
'''

import logging
import math
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import traceback

from Bio.PDB.PDBParser import PDBParser
from tqdm import tqdm

from SNDG import init_log, mkdir
from SNDG.Structure.ChainSplitter import ChainSplitter
from SNDG.Structure.FPocket import FPocket,fpocket_properties_map
from SNDG.Structure.PDBdb import *
from SNDG.Structure.PDBs import PDBs
from SNDG.Structure.QMean import QMean

init_log(rootloglevel=logging.INFO)
_log = logging.getLogger(__name__)

pocket_prop_map = {v:k for k,v in fpocket_properties_map.items()}



def process_chain(pdb_path, code, chain_id, pdb_model,props):
    chain_pdb_path = cs.make_pdb(pdb_path, code, chain_id, overwrite=True)
    qm = QMean.assesment(chain_pdb_path)
    residues_qm = qm["residues"]

    for k, v in residues_qm.items():
        chain, resid, resname = k.split("_")
        r = Residue.select().where(
            Residue.pdb == pdb_model & Residue.chain == chain & Residue.resid == int(resid)).first()
        for prop, val in v.items():
            prop_model = props["qr_" +  prop]
            if not math.isnan(val):
                ResidueProperty(residue=r, property=prop_model, value=val).save()

    del qm["residues"]
    for k, v in qm.items():
        prop_model = Property.select().where(Property.name == ("q_" + prop)).first()
        ChainProperty(pdb=pdb_model, chain=chain_id, name=k, value=v)#.save()
    res = FPocket(chain_pdb_path, chains_dir).hunt_pockets()
    for pocket in res.pockets:
        rs = ResidueSet(name="ChainPocket%i" % pocket.pocket_num, pdb=pdb_model)
        # rs.save()
        for k, v in pocket.properties.items():
            prop = "f_" + pocket_prop_map[k]
            prop_model = Property.select().where(Property.name == prop).first()
            ResidueSetProperty(residue_set=rs, property=prop_model, value=v)#.save()
    res.delete_dir()


if __name__ == '__main__':

    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-p", "--dbpass", required=True)
    parser.add_argument("-i", "--pdb_dir", default="/data/databases/pdb/")
    parser.add_argument("-db", "--dbname", default="pdbdb")
    parser.add_argument("-u", "--dbuser", default="root")

    args = parser.parse_args()
    from peewee import MySQLDatabase

    mysql_db = MySQLDatabase(args.dbname, user=args.dbuser, password=args.dbpass)

    sqldb.initialize(mysql_db)

    pdb_utils = PDBs(pdb_dir=args.pdb_dir)
    props = {x.name:x for x in Property.select()}
    pdbs = list(pdb_utils)
    with tqdm(pdbs) as pbar:
        for (code, pdb_path) in pbar:

            pdb_model = PDB.select().where(PDB.code == code).first()

            p = PDBParser(PERMISSIVE=True, QUIET=True)
            try:
                for chain in p.get_structure(code, pdb_path).get_chains():
                    chains_dir = args.pdb_dir + "/chains/" + code[1:3] + "/"
                    mkdir(chains_dir)
                    cs = ChainSplitter(chains_dir)
                    process_chain(pdb_path, code, chain.id, pdb_model,props)


            except Exception as ex:
                traceback.print_stack()
                _log.error(code + ": " + str(ex))
