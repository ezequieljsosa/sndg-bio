import os
import sys
import pandas as pd
from itertools import combinations
import Bio.SeqIO as bpio

class Primers:
    pass


if __name__ == "__main__":
    import argparse
    import Bio.SeqIO as bpio
    import os
    import fileinput
    import subprocess as sp
    import traceback

    parser = argparse.ArgumentParser(description='primer utils')
    # parser.add_argument('positional', action="store")
    # parser.add_argument('-p', '--primers', required=True)
    # parser.add_argument('-i', '--input_fasta', required=True)
    # parser.add_argument('-o', '--output', help="output_directory", default="./")
    # parser.add_argument('-o', '--output', help="output_directory", default="./")

    subparsers = parser.add_subparsers(help='commands', description='valid subcommands', required=True, dest='command')

    amplicon = subparsers.add_parser('amplicon', help='List contents')
    amplicon.add_argument('-i', '--blast_result', required=True)
    amplicon.add_argument('-s', '--sequences', required=True)
    # amplicon.add_argument('-o', '--output', help="output_directory", default="./")
    # subseq = subparsers.add_parser('hits', help='List contents')
    # subseq.add_argument('-o', '--output', help="output_directory", default="./")

    args = parser.parse_args()

    if args.command == "amplicon":
        df_blast = pd.read_table(args.blast_result,index_col=False,
                                 names=["query", "subject", "ident", "alnsize", "gap", "x", "qstart", "qend", "sstart",
                                        "send"],sep="\t")
        # df_blast["subject"] = [x.split("|")[-1] for x in df_blast["subject"]]
        seqs = bpio.to_dict( bpio.parse(args.sequences,"fasta"))

        for contig, df in df_blast.groupby("subject"):
            if len(df) > 1:

                for hit1,hit2 in combinations([x for _,x in df.iterrows()],2):
                    # print("----------------------------------------")
                    # print(hit1)
                    # print(hit2)
                    # print("----------------------------------------")
                    if hit1["query"] != hit2["query"]:
                        h1s = 1 if hit1.sstart < hit1.send else -1
                        h2s = 1 if hit2.sstart < hit2.send else -1
                        h1q = 1 if hit1.qstart < hit1.qend else -1
                        h2q = 1 if hit2.qstart < hit2.qend else -1

                        h1s = h1s * h1q
                        h2s = h2s * h2q
                        amplifica = False
                        reverso = False
                        if h1s != h2s:
                            if hit1.sstart < hit2.sstart:
                                if  (h1s == 1 and h2s == -1):
                                    amplifica = True
                            else:
                                reverso = True
                                if  (h1s == -1 and h2s == 1):
                                    amplifica = True
                        # print(f'{hit1["query"]} {hit2["query"]} {contig}')
                        # print(h1s,h2s,reverso,amplifica)
                        if amplifica:
                            min_pos = min([hit1.sstart,hit1.send,hit2.sstart,hit2.send] )
                            max_pos = max([hit1.sstart,hit1.send,hit2.sstart,hit2.send] )
                            print(f'{hit1["query"]} {hit2["query"]} {contig}-{min_pos}:{max_pos} reverso: {reverso}'  )
                            print( contig, len(seqs [contig][min_pos:max_pos].seq))
    # if args.primers == "-":
    #     for primer in fileinput.input("-"):
    #             pass
