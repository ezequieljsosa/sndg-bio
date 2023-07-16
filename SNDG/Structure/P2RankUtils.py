import numpy as np
import pandas as pd
from SNDG import Struct
import Bio.SeqIO as bpio
from Bio.PDB.PDBParser import PDBParser


def residues_near_drug(drug_centroid, aa_residues):
    residues_near = []
    for r in aa_residues:
        for a in list(r):
            dist = a - Struct(coord=drug_centroid)
            if dist > 20:
                break
            if dist < 10:
                residues_near.append(r)
                break
    return residues_near


def centeroid(arr):
    length = len(arr)
    sum_x = np.sum([x.coord[0] for x in arr])
    sum_y = np.sum([x.coord[1] for x in arr])
    sum_z = np.sum([x.coord[2] for x in arr])
    return sum_x / length, sum_y / length, sum_z / length


class P2RankUtils:

    @staticmethod
    def ligands_distance(pdb_name, model, contact_dist=6, exclude_drugs=[]):

        drug_molecules = []

        for res_obj in model.get_residues():
            if res_obj.id[0].strip():
                if res_obj.id[0] not in exclude_drugs:
                    drug_molecules.append(res_obj)

        for chain_obj in list(model):
            aa_residues = []
            chain = chain_obj.id
            for res_obj in chain_obj.get_residues():
                if not res_obj.id[0].strip():
                    aa_residues.append(res_obj)

            for res_drug_obj in drug_molecules:
                drug_centroid = centeroid(list(res_drug_obj))
                near_residues = residues_near_drug(drug_centroid, aa_residues)
                for drug_atom in list(res_drug_obj):
                    for near_residue in near_residues:
                        for residue_atom in list(near_residue):
                            distance = (residue_atom - drug_atom)
                            if distance > 20:
                                break
                            if distance < contact_dist:
                                fields = [pdb_name, chain, near_residue.id[1], near_residue.resname,
                                          residue_atom.serial_number,
                                          res_drug_obj.id[1], res_drug_obj.resname, drug_atom.serial_number, distance]
                                yield fields

    @staticmethod
    def parse_p2rank(prank_result_path):
        df = pd.read_csv(prank_result_path)
        df.columns = [x.strip() for x in df.columns]
        df["pocket"] = [x.strip() for x in df["name"]]
        df["residue_ids"] = [set([y.strip() for y in x.strip().split(" ")]) for x in df["residue_ids"]]
        df["surf_atom_ids"] = [set([y.strip() for y in x.strip().split(" ")]) for x in df["surf_atom_ids"]]
        return df

    @staticmethod
    def parse_dist(dist_result_path):
        df = pd.read_table(dist_result_path, index_col=False, names=["pdb_name", "chain", "res_id", "resname",
                                                                     "atom_serial_number",
                                                                     "ligand_id", "ligand_name", "ligand_atom_serial",
                                                                     "distance"])
        return df

    @staticmethod
    def intersect_pocket_ligand(df_p2rank, df_dist, solvent_cutoff=10, ligand_res_dist_cutoff=6):
        df_dist2 = df_dist[(df_dist.ligand_name != "HOH") & (df_dist.distance < ligand_res_dist_cutoff)]
        ligands_dict = df_dist2.groupby("ligand_name").agg({"ligand_id": lambda x: len(set(x))}).to_dict()["ligand_id"]
        exclude = [ligand_id for ligand_id, res_counts in ligands_dict.items() if res_counts > solvent_cutoff]
        df_dist2 = df_dist2[~df_dist2.ligand_name.isin(exclude)]
        ligand_residues = {}
        for ligand_id, df_ligand in df_dist2.groupby("ligand_id"):
            key = df_ligand.iloc[0].ligand_name + "_" + str(ligand_id)
            residues = [r.chain + "_" + str(r.res_id) for _, r in df_ligand.iterrows()]
            residues = set(residues)
            ligand_residues[key] = {"real": list(residues)}
            for _, p in df_p2rank[["pocket", "residue_ids"]].iterrows():

                intersect = list(p.residue_ids & residues)
                if intersect:
                    ligand_residues[key][str(p.pocket)] = intersect


        import json

        print(json.dumps(ligand_residues, indent=4))

    @staticmethod
    def pockets_to_ppsalignpocs(pdb1, prank_result1, pdb2, prank_result2):
        """
        https://seq2fun.dcmb.med.umich.edu//PPS-align/POC_format.pdf
        :return:
        """

        # with open() as h:

    @staticmethod
    def distance_intersect():
        pass


"""

"""

if __name__ == "__main__":
    # print(P2RankUtils.parse_p2rank("/tmp/pepe/103m/103m.pdb.gz_predictions.csv"))
    df_prank = P2RankUtils.parse_p2rank("/tmp/pepe/103m/103m.pdb.gz_predictions.csv")
    df_dist = P2RankUtils.parse_dist("/tmp/pepe/103m.dist")

    P2RankUtils.intersect_pocket_ligand(df_prank, df_dist)

    """
    import argparse

    parser = argparse.ArgumentParser(description='Utils for vcfs with multiple samples')

    subparsers = parser.add_subparsers(help='commands', description='valid subcommands', dest='command', required=True)
    cmd = subparsers.add_parser('res_near_lig', help='annotates a list of residues')
    cmd.add_argument('pdb_name', help='4 letter pdb code')
    cmd.add_argument('pdb_path', help='pdb_file_path')
    cmd.add_argument('--contact_dist', default=6, help='distance cutoff between a ligand and residue atom')
    cmd.add_argument('--include_hoh', action="store_true")

    cmd = subparsers.add_parser('phylo', help='prepare vcf for phylogeny')
    cmd.add_argument('vcf', help='joint vcf file')

    args = parser.parse_args()

    if args.command == "res_near_lig":
        p = PDBParser(PERMISSIVE=1, QUIET=1)
        model = list(p.get_structure(args.pdb_name, args.pdb_path))[0]
    """
