'''
Created on Nov 2, 2017

@author: eze
'''

import json
import os
import logging
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from collections import defaultdict
from tqdm import tqdm

import Bio.SeqIO as bpio
from Bio.PDB.PDBParser import PDBParser

from SNDG.Structure.PDBs import PDBs

from SNDG.Structure.ChainSplitter import ChainSplitter, SelectResidues
from SNDG.Structure.FPocket import FPocket
from SNDG.Structure.QMean import QMean




from SNDG import init_log,mkdir

init_log(rootloglevel=logging.INFO)
_log = logging.getLogger(__name__)

if __name__ == '__main__':

    parser = ArgumentParser( formatter_class=RawDescriptionHelpFormatter)


    parser.add_argument("-i", "--data_path", default='/data/databases/pdb/')
    parser.add_argument("-o", "--output_path", default='/data/databases/pdb/processed/domain_analisis')




    args = parser.parse_args()

    domains = defaultdict(lambda :[])
    for seq in bpio.parse(args.data_path + "/processed/domains.fasta","fasta"):
        domains["_".join(seq.id.split("_")[0:2]) ].append(seq.id.split("_"))
        
         
    for (code, pdb_path) in tqdm(PDBs()):

        p = PDBParser(PERMISSIVE=True, QUIET=True)
        try:
            for chain in p.get_structure(code, pdb_path).get_chains():
                chains_dir = args.output_path + "/chains/" + code[1:3] + "/"
                mkdir(chains_dir)
                cs = ChainSplitter(chains_dir)
                chain_pdb_path = cs.make_pdb(pdb_path,code,chain.id,overwrite=True)
                res = FPocket(chain_pdb_path, chains_dir).hunt_pockets()
                res.save(chains_dir + code + "_chain_" + chain.id + ".json")
                qm = QMean.analize(chain_pdb_path)
                with open(chain_pdb_path + ".qmean","w") as h:
                    json.dump(qm,h)
                for (_,_,res_start,res_end,dn,dn_start,dn_end) in domains[code + "_" + chain.id]:
                    #1r9d_A_2_787_PF02901.14_8_648     
                    try:
                        domains_dir = args.output_path + "/domains/" + code[1:3] + "/"
                        mkdir(domains_dir)
                        dn_start = int(dn_start)
                        dn_end   = int(dn_end)
                        cs.filter = SelectResidues(chain.id,{y:1 for y in [x.id[1] for x in chain.get_residues()][dn_start:dn_end]})
                        chain_pdb_path = cs.make_pdb(pdb_path,code,chain.id,overwrite=True)
                        res = FPocket(chain_pdb_path, domains_dir).hunt_pockets()
                        file_name = domains_dir + "/" + code + "_chain_" + \
                                    chain.id + "_".join([dn,str(dn_start),str(dn_end)])   +    ".json"
                        res.save(file_name)
                        res.delete_dir()
                        qm = QMean.analize(chain_pdb_path)
                        with open(file_name + ".qmean","w") as h:
                            json.dump(qm,h)
                    except Exception as ex:
                        _log.error(  "_".join([code,chain.id,res_start,res_end,dn,str(dn_start),str(dn_end)]) + ": " + str(ex))
        except Exception as ex:
            _log.error(  code + ": " + str(ex))
            
            
            
            
            