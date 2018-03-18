'''
Created on Nov 2, 2017

@author: eze
'''
from _collections import defaultdict

import Bio.SeqIO as bpio

from BiaStructure.IO import PDBsIterator
from Bio.PDB.PDBParser import PDBParser
from BiaStructure.IO.ChainSplitter import ChainSplitter, SelectResidues
from BiaStructure.Programs.fpocket import FPocket
from BiaStructure.Programs.QMean import QMean
import json
import os


if __name__ == '__main__':
    
    domains = defaultdict(lambda :[])
    for seq in bpio.parse("/data/databases/pdb/processed/domains.fasta","fasta"):         
        domains["_".join(seq.id.split("_")[0:2]) ].append(seq.id.split("_"))
        
         
    for (idx,(code, pdb_path)) in enumerate(PDBsIterator()):
        print str(idx) + " " + code        
        p = PDBParser(PERMISSIVE=True, QUIET=True)
        try:
            for chain in p.get_structure(code, pdb_path).get_chains():
                cs = ChainSplitter("/tmp/")
                chain_pdb_path = cs.make_pdb(pdb_path,code,chain.id,overwrite=True)
                res = FPocket(chain_pdb_path, "/tmp/").hunt_pockets()
                res.save("/tmp/" + code + "_chain_" + chain.id + ".json")
                
                for (_,_,res_start,res_end,dn,dn_start,dn_end) in domains[code + "_" + chain.id]:
                    #1r9d_A_2_787_PF02901.14_8_648     
                    try:    
                        dn_start = int(dn_start)
                        dn_end   = int(dn_end)
                        cs.filter = SelectResidues(chain.id,{y:1 for y in [x.id[1] for x in chain.get_residues()][dn_start:dn_end]})
                        chain_pdb_path = cs.make_pdb(pdb_path,code,chain.id,overwrite=True)
                        res = FPocket(chain_pdb_path, "/tmp/").hunt_pockets()
                        file_name = os.path.dirname(pdb_path) + "/" + code + "_chain_" + chain.id + "_".join([dn,str(dn_start),str(dn_end)])   +    ".json"
                        res.save(file_name)
                        res.delete_dir()
                        qm = QMean().analize(chain_pdb_path)
                        with open(file_name + ".qmean","w") as h:
                            json.dump(qm,h)
                    except Exception as ex:
                        print "error en " +  "_".join([code,chain.id,res_start,res_end,dn,str(dn_start),str(dn_end)]) + ": " + str(ex)
        except Exception as ex:
            print "error en " +  code + ": " + str(ex)
            
            
            
            
            