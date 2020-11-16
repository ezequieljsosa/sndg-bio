"""

"""

import os
from SNDG.WebServices import download_file
from SNDG.WebServices.Uniprot import Uniprot
import Bio.SeqIO as bpio
import gzip
from tqdm import  tqdm
class Reactome:
    DEFAULT_UNIP2REACTIONS = "https://reactome.org/download/current/UniProt2ReactomeReactions.txt"

    @staticmethod
    def download_fasta(url_unip2reactions=DEFAULT_UNIP2REACTIONS, outdir="/data/databases/reactome/", ovewrite=False):
        unip_utils = Uniprot()
        assert os.path.exists(outdir), f'{outdir} does not exists'
        reactome_map_file = outdir + "/UniProt2ReactomeReactions.txt"
        if ovewrite or not os.path.exists(reactome_map_file):
            download_file(url_unip2reactions, reactome_map_file, ovewrite=ovewrite)
        else:
            sys.stderr.write(f'{reactome_map_file} already exists')
        with open(reactome_map_file) as hr, gzip.open("seqs.fasta.gz", "wt") as hw:
            for line in tqdm(hr):
                if not line.startswith("#"):
                    unip, reactome, url_path, description = line.split("\t")[:4]
                    record = unip_utils.download_and_load_seqrecord(unip,format=".fasta")
                    record.name = ""
                    record.description = description + "||" + unip
                    record.id = reactome
                    bpio.write(record, hw, "fasta")


if __name__ == "__main__":
    import argparse
    import Bio.SeqIO as bpio
    import os
    import fileinput
    import subprocess as sp
    import traceback
    import sys

    parser = argparse.ArgumentParser(description='Reactome utils')
    subparsers = parser.add_subparsers(help='commands', description='valid subcommands', required=True, dest='command')

    download_cmd = subparsers.add_parser('download', help='download data from reactome')
    download_cmd.add_argument('--url', help="download url for uniprot to reactions map",
                              default=Reactome.DEFAULT_UNIP2REACTIONS)
    download_cmd.add_argument('-o', '--outdir', help="output_directory", default="./")
    download_cmd.add_argument('--force', action="store_true")
    download_cmd.add_argument('-f', '--format', default="fasta", choices=["fasta"], help="download format")

    args = parser.parse_args()
    if args.command == "download":
        sys.stderr.write(f'creating Reactome fasta...\n')
        Reactome.download_fasta(url_unip2reactions=args.url, outdir=args.outdir, ovewrite=args.force)
