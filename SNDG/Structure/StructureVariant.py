import sys

import Bio.SeqIO as bpio

from SNDG import execute, docker_wrap_command
from SNDG.Comparative.MSAMap import MSAMap
from SNDG.Structure.PDBs import PDBs, PDBFile
from SNDG.Structure.StructureAnnotator import StructureAnnotator
from collections import defaultdict


class StructureVariant:
    def __init__(self, pdb_dir="/data/databases/pdb/"):
        self.utils = PDBs(pdb_dir)
        self.seqs_path = "/tmp/seq.faa"
        self.aln_path = "/tmp/msa.faa"
        self.ref_seq = None
        self.pdbfile = None
        self.pdb_data = defaultdict(dict)

    def load_msa(self, input_sequence, pdb_code, pdb_chain=None):
        pdb_code = pdb_code.lower()
        self.utils.update_pdb(pdb_code)
        self.ref_seq = bpio.read(input_sequence, "fasta")
        self.pdbfile = PDBFile(pdb_code, self.utils.pdb_path(pdb_code))
        with open(self.seqs_path, "w") as h:
            bpio.write(self.ref_seq, h, "fasta")
            bpio.write(self.pdbfile.seq(selected_chain=pdb_chain), h, "fasta")

        cmd = docker_wrap_command(f'mafft --quiet --localpair --maxiterate 1000 {self.seqs_path} > {self.aln_path} ')
        execute(cmd)

        self.msa = MSAMap.from_msa(self.aln_path)
        self.res_map = self.pdbfile.residues_map(pdb_chain)

    def residues_from_pos(self, pos):
        pos_data = []
        for sample in self.msa.samples():
            if sample != self.ref_seq.id:
                pdb, chain = sample.split("_")[:2]
                if self.msa.exists_pos(self.ref_seq.id, pos, sample):
                    msa_pos = self.msa.pos_seq_msa_map[self.ref_seq.id][pos]
                    sample_pos = self.msa.pos_from_seq(self.ref_seq.id, pos, sample)
                    line = {
                        "pos": pos + 1,
                        "ref": self.msa.seqs[self.ref_seq.id][msa_pos],
                        "alt": self.msa.seqs[sample][msa_pos],
                        "pdb": pdb,
                        "chain": chain,
                        "resid": str(self.res_map[chain][sample_pos][1]),
                        "icode": str(self.res_map[chain][sample_pos][2]),
                        "pdb_pos": sample_pos
                    }
                    pos_data.append(line)
        return pos_data

    def residues_from_aln_seq(self, input_sequence, pdb_code, pdb_chain=None):
        self.load_msa(input_sequence, pdb_code, pdb_chain)
        variants = [(k, v) for k, v in
                    sorted(self.msa.variants(self.ref_seq.id).items(), key=lambda x: int(x[0].split("_")[1]))]

        output = []
        for ref_pos, alt_samples in variants:
            ref, pos = ref_pos.split("_")
            pos = int(pos)
            for alt, samples in alt_samples.items():
                if alt != self.msa.gap_code:
                    pos_data = self.residues_from_pos(pos)
                    output += pos_data
        return pd.DataFrame(output)

    def annotate_resid(self, pdb: str, resid: str, structure_annotator: StructureAnnotator):
        pdb = pdb.lower()
        data = {}
        if pdb not in self.pdb_data:
            self.load_pdb_ann(pdb, structure_annotator)

        if str(resid) in self.pdb_data[pdb]["binding"]:
            data["binding"] = self.pdb_data[pdb]["binding"][str(resid)]
        if str(resid) in self.pdb_data[pdb]["pockets"]:
            data["pockets"] = self.pdb_data[pdb]["pockets"][str(resid)]
        return data

    def load_pdb_ann(self, pdb, structure_annotator:StructureAnnotator):
        binding_data = structure_annotator.load_pdb_binding_data(pdb)
        binding_dict = defaultdict(list)
        for site in binding_data:
            for site_res in site["site_residues"]:
                res = str(site_res["residue_number"]) + (site_res.get("author_insertion_code", "") or "")
                binding_dict[res].append({
                    "site_id": site["site_id"],
                    "details": site["details"],
                    "ligands": [{c: x[c] for c in ["chain_id", "author_residue_number", "chem_comp_id"]} for x in
                                site["site_residues"] if
                                x["chem_comp_id"] in binding_dict and (x["chem_comp_id"] != "HOH")]
                })
        self.pdb_data[pdb]["binding"] = binding_dict
        pockets_data = structure_annotator.load_pdb_pocket(pdb,self.utils.pdb_dir)
        pockets_dict = defaultdict(list)
        for pocket in pockets_data:
            for residue in set(pocket["residues"]):
                pockets_dict[residue].append({"pocket_num": pocket["number"],
                                                               "druggabilitty": pocket["properties"][
                                                                   'Druggability Score']})
        self.pdb_data[pdb]["pockets"] = dict(pockets_dict)

    def annotate_residue_list(self, df, structure_annotator:StructureAnnotator):
        """
        
        :param df: columns=["pdb", "chain", "resid", "alt", "ref", "pos"] or generated by residues_from_aln_seq
        :return: 
        """

        output = {}
        for i, r in df.iterrows():
            output[f'{r.pdb}_{r.chain}_{r.resid}_{r.alt}'] = self.annotate_resid(r.pdb, str(r.resid),structure_annotator)

        return output


if __name__ == "__main__":
    import argparse
    import fileinput
    import pandas as pd
    import json

    parser = argparse.ArgumentParser(description='residue variants from sequence')

    subparsers = parser.add_subparsers(help='commands', description='valid subcommands', dest='command')
    seqcmd = subparsers.add_parser('residues', help='Get the residues map to a query sequence')

    seqcmd.add_argument('-i', '--input_sequence', default="-")
    seqcmd.add_argument('-db', '--pdb_data', default="/data/databases/pdb/")
    seqcmd.add_argument('-s', '--pdb_code', type=lambda s: s.lower())
    seqcmd.add_argument('-c', '--pdb_chain', type=lambda s: s.upper(), default=None)

    anncmd = subparsers.add_parser('ann', help='annotates a list of residues')
    anncmd.add_argument('-i', '--input_table', default="-", help='tsv with "pdb" "chain" "resid" to annotate')
    anncmd.add_argument('-db', '--pdb_data', default="/data/databases/pdb/")

    args = parser.parse_args()

    sv = StructureVariant(args.pdb_data)
    sa = StructureAnnotator()

    if args.command == "residues":
        df = sv.residues_from_aln_seq(fileinput.input(args.input_sequence), args.pdb_code, pdb_chain=args.pdb_chain)
        df.to_csv(sys.stdout, columns=["pdb", "chain", "resid", "alt", "ref", "pos", "pdb_pos"], 
                  index=False,  sep="\t")

    if args.command == "ann":
        from io import StringIO

        data = []
        for l in fileinput.input(args.input_table):
            data.append(l)
        df = pd.read_csv(StringIO("\n".join(data)), sep="\t")
        cols = set(df.columns)
        assert len(cols & set(["pdb", "chain", "resid"])) == 3
        for k,ann in sv.annotate_residue_list(df, sa).items():
            print(json.dumps({"residue":k,"ann":ann}, indent=4, sort_keys=True))
            print("###")
