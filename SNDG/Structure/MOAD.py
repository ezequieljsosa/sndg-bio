# https://bindingmoad.org/
from typing import BinaryIO
from collections import defaultdict


class MOAD:
    """
    https://bindingmoad.org/
    """

    @staticmethod
    def parse(moad_handle: BinaryIO):
        # family -> pdb -> compound -> residues
        result = {}
        comp_smiles_dict = {}
        for line in moad_handle:
            vec = line.split(",")
            if vec[0]:
                family = vec[0]
                result[family] = {}
            if vec[2]:
                pdb = vec[2]
                result[family][pdb] = []
            if vec[3]:
                hetatom = vec[3]
                compound, chain, resid = hetatom.split(":")
                status = vec[4]
                comp_data = {"status": status, "compound": compound, "chain": chain, "resid": resid}
                result[family][pdb].append(comp_data)
                if vec[5].startswith("K"):
                    comp_data["affinity"] = {"metric": vec[5], "relation": vec[6], "measure": vec[7], "units": vec[8]}
                if vec[-2]:
                    comp_smiles_dict[compound] = vec[-2]
        return result, comp_smiles_dict


if __name__ == "__main__":
    import argparse
    import json
    import sys

    parser = argparse.ArgumentParser(description='create json from moad csv')
    parser.add_argument('moad',  help='moad path')
    parser.add_argument('--smiles_out',  help='smiles output file',default=None)
    args = parser.parse_args()
    with open(args.moad) as h:
        result, comp_smiles_dict = MOAD.parse(h)
        json.dump(result, sys.stdout, indent=4, sort_keys=True)
    if args.smiles_out:
        with open(args.smiles_out) as h:
            json.dump(comp_smiles_dict, args.smiles_out, indent=4, sort_keys=True)

