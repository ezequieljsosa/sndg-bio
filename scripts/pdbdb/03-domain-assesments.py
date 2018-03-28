'''
Created on Nov 2, 2017

@author: eze
'''

import json
import os
import logging
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from collections import defaultdict

import math
from peewee import MySQLDatabase
from tqdm import tqdm

import Bio.SeqIO as bpio
from Bio.PDB.PDBParser import PDBParser

from SNDG.Structure.PDBs import PDBs

from SNDG.Structure.ChainSplitter import ChainSplitter, SelectResidues
from SNDG.Structure.FPocket import FPocket
from SNDG.Structure.QMean import QMean

from SNDG.Structure.PDBdb import *

from SNDG import init_log, mkdir

init_log(rootloglevel=logging.INFO)
_log = logging.getLogger(__name__)

mysql_db = MySQLDatabase('pdbdb', user="root", password="mito")

sqldb.initialize(mysql_db)


def process_chain(pdb_path, code, chain_id, pdb_model):
    chain_pdb_path = cs.make_pdb(pdb_path, code, chain_id, overwrite=True)
    qm = QMean.assesment(chain_pdb_path)
    residues_qm = qm["residues"]

    for k,v in residues_qm.items():
        chain, resid, resname = k.split("_")



        r = Residue(pdb=pdb_model,chain=chain,resname=resname,resid=int(resid))
        for prop,val in v.items():
            if not math.isnan(val):
                ResidueProperty(residue=r,name=prop,value=val).save()

    del qm["residues"]
    for k, v in qm.items():
        ChainProperty(pdb=pdb_model, chain=chain.id, name=k, value=v).save()
    res = FPocket(chain_pdb_path, chains_dir).hunt_pockets()
    for pocket in res.pockets:
        rs = ResidueSet(name="ChainPocket%i" % pocket.pocket_num, pdb=pdb_model)
        rs.save()
        for k, v in pocket.properties.items():
            ResidueSetProperty(residue_set=rs, name=k, value=v).save()
    res.delete_dir()


def process_domain(domains_dir, chain, dn_start, dn_end, pdb_model):
    mkdir(domains_dir)

    cs.filter = SelectResidues(chain.id, {y: 1 for y in
                                          [x.id[1] for x in chain.get_residues()][dn_start:dn_end]})
    domain_pdb_path = cs.make_pdb(pdb_path, code, chain.id, overwrite=True)
    res = FPocket(domain_pdb_path, domains_dir).hunt_pockets()
    for pocket in res.pockets:
        rs = ResidueSet(name="DomainPocket%i" % pocket.pocket_num, pdb=pdb_model)
        rs.save()
        for k, v in pocket.properties.items():
            ResidueSetProperty(residue_set=rs, name=k, value=v).save()
    res.delete_dir()

    qm = QMean.assesment(domain_pdb_path)
    residues_qm = qm["residues"]
    del qm["residues"]
    for k, v in qm.items():
        ChainProperty(pdb=pdb_model, chain=chain.id, name=k, value=v).save()


if __name__ == '__main__':

    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--data_path", default='/data/databases/pdb/')
    parser.add_argument("-o", "--output_path", default='/data/databases/pdb/processed/domain_analisis')

    args = parser.parse_args()

    domains = defaultdict(lambda: [])
    for seq in bpio.parse(args.data_path + "/processed/domains.fasta", "fasta"):
        domains["_".join(seq.id.split("_")[0:2])].append(seq.id.split("_"))

    for (code, pdb_path) in tqdm(PDBs(pdb_dir=args.data_path)):

        pdb_model = PDB(code=code)
        pdb_model.save()

        p = PDBParser(PERMISSIVE=True, QUIET=True)
        try:
            for chain in p.get_structure(code, pdb_path).get_chains():
                chains_dir = args.output_path + "/chains/" + code[1:3] + "/"
                mkdir(chains_dir)
                cs = ChainSplitter(chains_dir)
                process_chain(pdb_path, code, chain.id, pdb_model)

                for (_, _, res_start, res_end, dn, dn_start, dn_end) in domains[code + "_" + chain.id]:
                    # 1r9d_A_2_787_PF02901.14_8_648
                    try:
                        domains_dir = args.output_path + "/domains/" + code[1:3] + "/"
                        dn_start = int(dn_start)
                        dn_end = int(dn_end)
                        process_domain(domains_dir, chain, dn_start, dn_end, pdb_model)
                    except Exception as ex:
                        _log.error(
                            "_".join([code, chain.id, res_start, res_end, dn, str(dn_start), str(dn_end)]) + ": " + str(
                                ex))
        except Exception as ex:
            _log.error(code + ": " + str(ex))
