from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import seq1
import Bio.SeqIO as bpio
from SNDG import mkdir

class PDBFile:

    def __init__(self, code, in_pdb):
        self.code = code
        self.struct = PDBParser(PERMISSIVE=1, QUIET=1).get_structure(self.code, in_pdb)

    def residues_map(self, selected_chain=None, standard_aa=True):
        rmap = {}
        for chain in self.struct.get_chains():
            if (not selected_chain) or (selected_chain == chain.id):
                residues = [x for x in chain.get_residues() if is_aa(x, standard=standard_aa)]
                rmap[chain.id] = {i: x.id[1] for i, x in enumerate(residues)}
        return rmap

    def seq(self, selected_chain=None, standard_aa=True):
        records = []
        for chain in self.struct.get_chains():
            if (not selected_chain) or (selected_chain == chain.id):
                residues = [x for x in chain.get_residues() if is_aa(x, standard=standard_aa)]
                if residues:
                    seq = "".join([seq1(x.resname) for x in residues])
                    start = str(residues[0].id[1])
                    end = str(residues[-1].id[1])
                    record = SeqRecord(id="_".join([self.code, chain.id, start, end]), description="", seq=Seq(seq))
                    records.append(record)
        return records


if __name__ == '__main__':
    import sys
    from SNDG import init_log
    import json
    import os
    import gzip

    init_log()

    import argparse
    from SNDG.Structure.PDBs import PDBs

    parser = argparse.ArgumentParser(description='PDB utils')

    parser.add_argument('--pdbdb', default=os.environ.get("PDBDB", None))

    subparsers = parser.add_subparsers(help='commands', description='valid subcommands', dest='command', required=True)

    cmd = subparsers.add_parser('seq', help='gets the sequence and residue map from a PDB code')
    cmd.add_argument('pdb', help="pdb_code")
    cmd.add_argument('--pdb_file', help="pdb_file", default=None)
    cmd.add_argument('--output_dir', help="output dir", default=None)

    args = parser.parse_args()

    if not args.pdbdb:
        parser.error(
            f'pdbdb parameter not defined. Either call the script with --pdbdb or  set the PDBDB environment variable ')

    pdbs = PDBs(pdb_dir=args.pdbdb)
    args.pdb = args.pdb.lower()

    if args.output_dir:
        pdbs_seq_path = f"{args.output_dir}/{args.pdb}.fasta.gz"
        pdbs_seqmap_path = f"{args.output_dir}/{args.pdb}_resmap.json.gz"

    else:

        pdbs_seq_path = pdbs.pdbseq_path(args.pdb)
        pdbs_seqmap_path = pdbs.pdbseqres_map_path(args.pdb)
        mkdir(os.path.dirname(pdbs_seq_path))

    if args.command == "seq":

        if args.pdb_file:
            pdb_path = args.pdb_file
        else:
            assert len(args.pdb) == 4, "PDB shoud be 4 letters long"
            pdb_path = pdbs.pdb_path(args.pdb.lower())
        assert os.path.exists(pdb_path), f'"{pdb_path}" does not exists...'

        pdb = PDBFile(args.pdb, pdb_path)
        with gzip.open(pdbs_seq_path, "wt") as hseq, gzip.open(
                pdbs_seqmap_path, "wt") as hmap:
            bpio.write(pdb.seq(), hseq, "fasta")
            json.dump(pdb.residues_map(), hmap)


        sys.exit(0)
