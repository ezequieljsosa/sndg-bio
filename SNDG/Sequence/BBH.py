"""

"""
import logging
from tqdm import tqdm
from collections import defaultdict
import Bio.SearchIO as bpsio
from SNDG.Sequence import read_blast_table
from SNDG import execute, docker_wrap_command

_log = logging.getLogger(__name__)


class BBH:

    @staticmethod
    def bbhs_from_blast(path_file1, path_file2, ident_threshold=80):
        """
        :param path_file1: blast result table format (6)
        :param path_file2: blast result table format (6)
        :return: list of bbh tuples
        """
        query_dict = defaultdict(lambda: {})
        for _, r in read_blast_table(path_file1).iterrows():

            if r.identity > ident_threshold:
                query_dict[r["query"]][r["hit"]] = 1

        bbhs = []
        for _, r in read_blast_table(path_file2).iterrows():
            if r.identity > ident_threshold:
                if r["query"] in query_dict[r["hit"]]:
                    bbhs.append((r["hit"], r["query"]))
        return bbhs

    @staticmethod
    def bbhs_from_lists(path_file1, path_file2, strict=False):
        """
        :param path_file1: blast result table format (6)
        :param path_file2: blast result table format (6)
        :return: list of bbh tuples
        """
        query_dict = defaultdict(lambda: {})
        for line in open(path_file1).readlines():
            try:
                x, y = line.strip().split("\t")[0:2]
                query_dict[x][y] = 1
            except:
                if strict:
                    raise Exception("parse Error: " + line)
                else:
                    _log.warn("parse Error: " + line)

        bbhs = []
        for line in open(path_file2).readlines():
            x, y = line.strip().split("\t")
            if x in query_dict[y]:
                bbhs.append([y, x])
        return bbhs


if __name__ == '__main__':
    import os
    import argparse
    import sys

    parser = argparse.ArgumentParser(description='Process BBH between 2 proteomes')
    parser.add_argument('-f1', '--fasta_file1', default=None)
    parser.add_argument('-f2', '--fasta_file2', default=None)

    # parser.add_argument('-d','--datadir', default=None)

    parser.add_argument('-b1', '--blast_result1', default=None)
    parser.add_argument('-b2', '--blast_result2', default=None)

    parser.add_argument('-i', '--ident_threshold', default=0.7)
    parser.add_argument('-e', '--evalue', type=float, default=1e-6)
    parser.add_argument('-qcov', '--query_cover', type=int, default=80)
    parser.add_argument('-scov', '--subject_cover', type=int, default=80)

    parser.add_argument('--force', action="store_true")

    # parser.add_argument('-outfmt','--ident_threshold', default=0.9)

    args = parser.parse_args()

    if args.fasta_file1 and args.fasta_file2:
        assert os.path.exists(args.fasta_file1), f'{args.fasta_file1} does not exists'
        assert os.path.exists(args.fasta_file2), f'{args.fasta_file2} does not exists'
    elif args.blast_result1 and args.blast_result2:
        assert os.path.exists(args.blast_result1), f'{args.blast_result1} does not exists'
        assert os.path.exists(args.blast_result2), f'{args.blast_result2} does not exists'
    else:
        sys.stderr.write(
            f'you must specify either fasta_file1 and fasta_file2 or blast_result1 and blast_result2 parameters\n')
        parser.print_help(sys.stderr)
        sys.exit(1)

    for force, blast_result, fasta_file1, fasta_file2 in [
        [args.force, args.blast_result1, args.fasta_file1, args.fasta_file2],
        [args.force, args.blast_result2, args.fasta_file2, args.fasta_file1]]:
        if not blast_result:
            blast_result = fasta_file1 + ".tbl"
        if force or not blast_result or not os.path.exists(blast_result):
            if force or not os.path.exists(fasta_file1 + ".dmnd"):
                cmd = f"diamond-aligner makedb  -d {os.path.abspath(fasta_file1)}  --in {os.path.abspath(fasta_file1)}"
                execute(docker_wrap_command(cmd))

            cmd = f"""diamond-aligner blastp -d {os.path.abspath(fasta_file1)}  -q {os.path.abspath(fasta_file2)} -e {args.evalue} \
                        --query-cover {args.query_cover} --subject-cover {args.subject_cover} """
            with open(blast_result, "w") as h:
                execute(docker_wrap_command(cmd), stdout=h)

    if not args.blast_result1:
        args.blast_result1 = args.fasta_file1 + ".tbl"
    if not args.blast_result2:
        args.blast_result2 = args.fasta_file2 + ".tbl"

    data = BBH.bbhs_from_blast(args.blast_result1, args.blast_result2, ident_threshold=args.ident_threshold)

    for x, y in data:
        print(x + "\t" + y)
