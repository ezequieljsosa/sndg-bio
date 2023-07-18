import pandas as pd
from tqdm import tqdm
import Bio.SeqIO as bpio
from Bio import AlignIO
from SNDG.Comparative.MSAMap import MSAMap
from SNDG.Structure.PDBs import PDBs

if __name__ == '__main__':
    import argparse
    import sys
    import os
    import subprocess as sp
    import json
    import tempfile
    import gzip
    from collections import defaultdict


    def is_valid_file(parser, arg):
        if not os.path.exists(arg):
            parser.error("The file %s does not exist!" % arg)
        else:
            return arg  # return an open file handle


    parser = argparse.ArgumentParser(description='annotates')

    parser.add_argument('--output_dir', nargs='?', default="./msa_results")
    parser.add_argument('--pdbdb', default=os.environ.get("PDBDB", None))
    parser.add_argument('--threads', default=1, type=int)

    subparsers = parser.add_subparsers(help='commands', description='valid subcommands', dest='command', required=True)
    cmd = subparsers.add_parser('addpdb', help='fix genebank file to be imported')
    cmd.add_argument('ref_accession', help="accession in the MSA to be used as a reference")
    cmd.add_argument('msa_file', help="MSA in fasta format", type=lambda x: is_valid_file(parser, x), metavar="FILE")
    cmd.add_argument('pdbs', nargs="+", help="MSA in fasta format")

    # cmd.add_argument('af_pdb_file', help="AlphaFold PDB file", type=lambda x: is_valid_file(parser, x), metavar="FILE")
    # cmd.add_argument('af_scores_file', help="AlphaFold PDB scores", type=lambda x: is_valid_file(parser, x),
    #                 metavar="FILE")
    cmd.add_argument('--pdbs_fasta', required=False, help="fasta file with the pdbs chain sequences")
    cmd.add_argument('--pdbs_resmap', required=False, help="JSON files with the mapping between the fasta " +
                                                           "sequence index and the residue ID")

    # cmd.add_argument('--distances', required=False, type=lambda x: is_valid_file(parser, x), metavar="FILE")

    cmd = subparsers.add_parser('addcsa', help='fix genebank file to be imported')
    cmd.add_argument('--csa', default="./csa_curated_data.csv", type=lambda x: is_valid_file(parser, x), metavar="FILE")
    cmd.add_argument('table', default=None, metavar="FILE", help="table with pdb,chain and resid columns")

    cmd = subparsers.add_parser('adddist', help='fix genebank file to be imported')
    cmd.add_argument('table', default=None, metavar="FILE", help="table with pdb,chain and resid columns")
    cmd.add_argument('--distfile', required=False, default=None, metavar="FILE",
                     help="table with pdb distances, only if 'pdb' is used")
    cmd.add_argument('--pdb', required=False, default=None, help="pdb code, only if 'distfile' is used ")
    cmd.add_argument('--filter', default=["HOH"], help="PDB ligand residues to be filtered")
    cmd.add_argument('--maxdist', default=6, help="max distance from the ligand")

    cmd = subparsers.add_parser('addtable')
    cmd.add_argument('basetable', type=lambda x: is_valid_file(parser, x), metavar="FILE")
    cmd.add_argument('addtable',  type=lambda x: is_valid_file(parser, x), metavar="FILE")
    cmd.add_argument('base_join_field')
    cmd.add_argument('add_join_field')
    cmd.add_argument('colum2add')

    args = parser.parse_args()

    if not args.pdbdb:
        parser.error(
            f'pdbdb parameter not defined. Either call the script with --pdbdb or  set the PDBDB environment variable ')

    pdbs = PDBs(pdb_dir=args.pdbdb)

    if args.command == "addtable":
        df1 = pd.read_table(args.basetable, sep="\t", index_col=False)
        df2 = pd.read_table(args.addtable, sep="\t", index_col=False)

        dfjoined = df1.merge(df2[ [args.add_join_field] + args.colum2add.split(",")],
                            left_on=args.base_join_field, right_on=args.add_join_field,
                             suffixes=('', '_second'), how="left")
        dfjoined.drop(columns=[(x + "_second") if (x + "_second") in dfjoined.columns else x
                               for x in args.add_join_field.split(",")],inplace=True)
        dfjoined.to_csv(sys.stdout, sep="\t", index=False)

    if args.command == "adddist":
        df = pd.read_table(args.table, sep="\t", index_col=False)
        df["pdb_chain_res"] = ["_".join([str(r.pdb), str(r.chain), str(r.resid)]) for _, r in
                               df[["pdb", "chain", "resid"]].iterrows()]

        if args.distfile:
            data = [[args.pdb, args.distfile]]
        else:
            pdb_list = set(df.pdb)
            data = [[pdb, pdbs.pdbdist_path(pdb)] for pdb in pdb_list]

        dist_dict = defaultdict(list)
        for pdb, distfile in data:
            dfdist = pd.read_table(distfile, compression='gzip', names=PDBs.DIST_COLS, sep="\t")
            dfdist = dfdist[~dfdist.ligname.isin(args.filter)][dfdist.dist < args.maxdist]
            dfdist["pdb_chain_res"] = ["_".join([r.pdb, r.chain, str(r.resid)]) for
                                       _, r in dfdist["pdb chain resid".split()].iterrows()]

            # "pdb chain resid resname atom ligid ligname atomlig dist"
            #  3ff5	B	23	GLU	499	106	HOH	1016	4.3496504

            for pdb_chain_res, df2 in dfdist.groupby("pdb_chain_res"):
                for ligid, df3 in df2.groupby("ligid"):
                    ligname = df3.iloc[0].ligname
                    pdb = df3.iloc[0].pdb
                    dist = df3["dist"].mean()
                    dist_dict[pdb_chain_res].append(f'{ligid}|{ligname}|{dist}')

        dist_dict = {k: " ".join(v) for k, v in dist_dict.items()}
        dfligs = pd.DataFrame(dist_dict.items(), columns=["pdb_chain_res", "ligs"])
        # dfdist["ligs"] = [dist_dict[x] for x in dfdist["pdb_chain_res"]]

        dfjoined = df.merge(dfligs,
                            left_on='pdb_chain_res', right_on='pdb_chain_res', suffixes=('', ''), how="left")

        dfjoined.drop(columns=["pdb_chain_res"], inplace=True)
        dfjoined.to_csv(sys.stdout, sep="\t", index=False)

    if args.command == "addcsa":

        df = pd.read_table(args.table, sep="\t", index_col=False)
        df["pdb_chain_res"] = ["_".join([str(r.pdb), str(r.chain), str(r.resid)]) for _, r in
                               df[["pdb", "chain", "resid"]].iterrows()]
        dfcsa = pd.read_csv(args.csa, on_bad_lines='skip')
        dfcsa["pdb_chain_res"] = ["_".join([str(r.PDB), str(r["chain/kegg compound"]), str(r["resid/chebi id"])]) for
                                  _, r in
                                  dfcsa[["PDB", "chain/kegg compound", "resid/chebi id"]].iterrows()]
        dfcsa["pdb_chain_res"] = dfcsa.pdb_chain_res.astype(str)

        dfcsa.rename(columns={"M-CSA ID": "csaid", "residue/reactant/product/cofactor": "ligand",
                              "role group": "role_group"}, inplace=True)

        df["pdb_chain_res"] = df.pdb_chain_res.astype(str)


        # pdb\tchain\tresid\tpdb_seq_pos\tpdb_aa\tref_pos\tmsa_aa\tmsa_pos
        # 'M-CSA ID', 'Uniprot IDs', 'PDB', 'EC',     'residue/reactant/product/cofactor', 'PDB code', 'chain/kegg compound',
        #        'resid/chebi id', 'function location/name', 'role', 'role type',     'role group'

        def valueset(seriesx):
            values = list(set(seriesx.dropna()))
            if len(values) == 1:
                return values[0]
            return values


        def setconcat(seriesx):
            return "|".join([str(x) for x in list(set(seriesx.dropna()))])


        # "ligand": valueset,
        dfcsa = dfcsa.groupby("pdb_chain_res").agg(
            {"csaid": "first", "role": setconcat, "role_group": setconcat}).reset_index()

        dfjoined = df.merge(dfcsa,
                            left_on='pdb_chain_res', right_on='pdb_chain_res', suffixes=('', '_csa'), how="left")

        dfjoined.drop(columns=["pdb_chain_res"], inplace=True)

        dfjoined.to_csv(sys.stdout, sep="\t", index=False)

    if args.command == "addpdb":

        if not args.pdbs_fasta:
            args.pdbs_fasta = tempfile.mkstemp()[1]
            resmap = {}

            with open(args.pdbs_fasta, "w") as htmpfasta:
                for pdb in args.pdbs:
                    pdb = pdb.lower()
                    with gzip.open(pdbs.pdbseq_path(pdb), "rt") as hseq1, gzip.open(pdbs.pdbseqres_map_path(pdb),
                                                                                    "rt") as hseq2:
                        htmpfasta.write(hseq1.read())
                        resmap[pdb] = json.load(hseq2)
        else:
            with open(args.pdbs_resmap) as h:
                resmap = json.load(h)

        alignments = bpio.to_dict(AlignIO.read(args.msa_file, "fasta"))

        if args.ref_accession not in alignments:
            sys.stderr.write(f'"{args.ref_accession}" is not in "{args.msa_file}"')
            sys.exit(1)
        if not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir)
            if not os.path.exists(args.output_dir):
                parser.error(f"could not create output dir: '{args.output_dir}'")

        cmd = f'mafft --thread {args.threads} --add {args.pdbs_fasta} --localpair {args.msa_file} > {args.output_dir}/msa_with_pdb.fasta'
        if not os.path.exists(f'{args.output_dir}/msa_with_pdb.fasta'):
            sp.call(cmd, shell=True)
        else:
            sys.stderr.write(
                f"WARNNIG '{args.output_dir}/msa_with_pdb.fasta' already exists, it will not be re-calculated\n")

        msa_map = MSAMap.from_msa(f'{args.output_dir}/msa_with_pdb.fasta')
        msa_map.init()

        with open(f'{args.output_dir}/msa_with_pdb.tbl', "w") as hw:
            hw.write(f'pdb\tchain\tresid\tpdb_seq_pos\tpdb_aa\tref_pos\tmsa_aa\tmsa_pos\n')

            for record in bpio.parse(args.pdbs_fasta, "fasta"):
                pdb, chain, resstart, resend = record.id.split("_")
                # print(record.id)
                chain_resmap = resmap[pdb][chain]
                for pdb_pos, pdb_aa in enumerate(str(record.seq)):
                    res_id = chain_resmap[str(pdb_pos)]
                    # print(pdb_pos,pdb_aa,res_id)
                    msa_pos = msa_map.pos_seq_msa_map[record.id][pdb_pos]
                    try:
                        ref_pos = msa_map.pos_from_seq(record.id, pdb_pos, args.ref_accession)
                    except ValueError:
                        ref_pos = "-"

                    ref_aa = msa_map.seqs[args.ref_accession][msa_pos]
                    hw.write(f'{pdb}\t{chain}\t{res_id}\t{pdb_pos}\t{pdb_aa}\t{ref_pos}\t{ref_aa}\t{msa_pos}\n')
