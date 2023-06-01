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
    def bbhs(query_best_hit_1, query_best_hit_2):
        """
        :param query_best_hit_1: dict of queries with best hit
        :param query_best_hit_2: dict of queries with best hit
        :return: bidirectional hit between inputs as a list of tuples -> (query, hit)
        """
        query_best_hit_1 = list(query_best_hit_1.items())
        query_best_hit_2 = list(query_best_hit_2.items())
        query_best_hit_2 = [t[::-1] for t in query_best_hit_2]
        bbhs = []
        for q in query_best_hit_1:
            if q in query_best_hit_2:
                bbhs.append((q))
        return bbhs        
    
    @staticmethod
    def best_hit(query_hits_dict):
        """
        :param query_hits_dict: dict of queries as key and hits as values including relevant blast values
        :return: dict of queries as key and only one hit as value (Best hit) -> {query: hit}
        """
        for query in query_hits_dict:
            max_identity = 0
            max_bitscore = 0
            min_evalue = 1
            selected_hit = ""
            for hit in query_hits_dict[query]:
                identity = query_hits_dict[query][hit][0]
                bitscore = query_hits_dict[query][hit][1]
                evalue = query_hits_dict[query][hit][2]
                if identity>max_identity and bitscore>max_bitscore and evalue<min_evalue:
                    max_identity = identity
                    max_bitscore = bitscore
                    min_evalue = evalue
                    selected_hit = hit
            query_hits_dict[query]=selected_hit
        return query_hits_dict

    @staticmethod
    def parse_blast(path_file, ident_threshold=80):
        """
        :param path_file: blast result table format (6)
        :return: dict of queries as key and hits as values including relevant blast values -> {query: {hit: (identity, bitscore, evalue), hit: (identity, bitscore, evalue)}}
        """
        query_hits_dict = defaultdict(lambda: {})
        for _, r in read_blast_table(path_file).iterrows():
            if r.identity > ident_threshold:
                query_hits_dict[r["query"]][r["hit"]] = (r.identity, r.bitscore, r.evalue)
        return query_hits_dict


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

    parser.add_argument('-i', '--ident_threshold', default=70)
    parser.add_argument('-e', '--evalue', type=float, default=1e-6)
    parser.add_argument('-qcov', '--query_cover', type=int, default=80)
    parser.add_argument('-scov', '--subject_cover', type=int, default=80)

    parser.add_argument('--force', action="store_true")

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

    query_hits_dict_1 = BBH.parse_blast(args.blast_result1, ident_threshold=args.ident_threshold)
    query_hits_dict_2 = BBH.parse_blast(args.blast_result2, ident_threshold=args.ident_threshold)
       
    query_best_hit_1 = BBH.best_hit(query_hits_dict_1)
    query_best_hit_2 = BBH.best_hit(query_hits_dict_2)

    bbhs = BBH.bbhs(query_best_hit_1, query_best_hit_2)

    for x, y in bbhs:
        print(x + "\t" + y)
        