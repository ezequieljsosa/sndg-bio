import os
import sys


class Primers:
    pass


if __name__ == "__main__":
    import argparse
    import Bio.SeqIO as bpio
    import os
    import fileinput
    import subprocess as sp
    import traceback

    parser = argparse.ArgumentParser(description='Download NCBI assemblies')
    # parser.add_argument('positional', action="store")
    # parser.add_argument('-p', '--primers', required=True)
    # parser.add_argument('-i', '--input_fasta', required=True)
    # parser.add_argument('-o', '--output', help="output_directory", default="./")
    # parser.add_argument('-o', '--output', help="output_directory", default="./")

    subparsers = parser.add_subparsers(help='commands', description='valid subcommands', required=True, dest='command')

    amplicon = subparsers.add_parser('amplicon', help='List contents')
    amplicon.add_argument('-o', '--output', help="output_directory", default="./")
    subseq = subparsers.add_parser('hits', help='List contents')
    subseq.add_argument('-o', '--output', help="output_directory", default="./")

    args = parser.parse_args()
    print(args)

    # if args.primers == "-":
    #     for primer in fileinput.input("-"):
    #             pass
