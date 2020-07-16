import os
import sys
import json
from collections import defaultdict

from SNDG import mkdir
from SNDG.Structure.FPocket import FPocket
from SNDG.Structure.PDBs import PDBs
from SNDG.WebServices.PDBsWS import PDBsWS


class StructureAnnotator:
    def __init__(self):
        self.pdb_data = defaultdict(dict)

    def load_pdb_binding_data(self, pdb):
        self.pdb_data[pdb]["ligands"] = PDBsWS.ligands(pdb)
        self.pdb_data[pdb]["binding"] = PDBsWS.binding(pdb)
        return self.pdb_data[pdb]["binding"]

    def load_pdb_pocket(self, pdb, pdb_dir="/data/databases/pdb/"):
        utils = PDBs(pdb_dir)
        if not os.path.exists(utils.pdb_pockets_path(pdb)):
            utils.update_pdb(pdb)
            fpocket = FPocket(utils.pdb_path(pdb))
            result = fpocket.hunt_pockets()
            mkdir(os.path.dirname(utils.pdb_pockets_path(pdb)))
            result.save(utils.pdb_pockets_path(pdb))
        with open(utils.pdb_pockets_path(pdb)) as h:
            result = json.load(h)

        self.pdb_data[pdb]["pockets"] = result
        return self.pdb_data[pdb]["pockets"]

    def load_pdb_ann_data(self, pdb):
        if pdb not in self.pdb_data:
            self.pdb_data[pdb] = {}
            self.load_pdb_binding_data(pdb)


if __name__ == "__main__":
    import argparse
    import fileinput

    parser = argparse.ArgumentParser(description='residue variants from sequence')

    parser.add_argument('-i', '--pdb_code', required=True, help='pdb code or list')
    parser.add_argument('-db', '--pdb_data', default="/data/databases/pdb/")
    parser.add_argument('-s', '--steps', nargs='*', default=["binding", "pocket"], choices=["binding", "pocket"])

    args = parser.parse_args()

    sv = StructureAnnotator()
    if os.path.exists(args.pdb_code) or (args.pdb_code == "-"):
        pdbs = fileinput.input(args.pdb_code)
    else:
        pdbs = args.pdb_code.split(",")
    for pdb in pdbs:
        pdb = pdb.lower().strip()
        if len(pdb) != 4:
            sys.stderr.write(f'{pdb} invalid pdb code')
            continue

        for step in args.steps:
            if step == "pocket":
                sv.load_pdb_pocket(pdb, args.pdb_data)
            if step == "binding":
                sv.load_pdb_binding_data(pdb)

        print({"pdb": pdb, "ann": sv.pdb_data[pdb]})

# https://www.ebi.ac.uk/pdbe/api/doc/pdb.html
# https://www.ebi.ac.uk/pdbe/api/mappings/:accession
# https://www.ebi.ac.uk/pdbe/api/pdb/entry/secondary_structure/:pdbid
# https://www.ebi.ac.uk/pdbe/api/pdb/entry/modified_AA_or_NA/:pdbid
# https://www.ebi.ac.uk/pdbe/api/pdb/entry/polymer_coverage/:pdbid/chain/:chainid
