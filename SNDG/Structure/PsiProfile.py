#!/usr/bin/env python3
"""

"""

import glob
import os
from tempfile import mkdtemp,mktemp
from xml.etree.ElementTree import ParseError

import Bio.SearchIO as bpsio
from SNDG import execute

import logging

_log = logging.getLogger(__name__)


class PsiProfile:
    search_result_fields = ["query", "qstart", "qend", "hit", "hstart", "hend", "evalue", "identity", "qaln", "haln"]

    @staticmethod
    def build_profile(seq_fasta, database, iterations, pssm_file, cpu, evalue=0.0001):
        cmd = f"psiblast -query {seq_fasta} -db {database} -num_threads {cpu} -out_pssm {pssm_file} -evalue {evalue} -num_iterations {iterations} 1>&2 2>/dev/null"
        execute(cmd)

    @staticmethod
    def profile_search(database, pssm_file, search_result, cpu=1, evalue=0.00001):
        cmd = f"psiblast -db {database} -in_pssm {pssm_file} -num_threads {cpu} -evalue {evalue}  -outfmt 5 -out {search_result} 1>&2"
        execute(cmd)
        try:
            search_result = list(bpsio.parse(search_result, "blast-xml"))
        except ParseError:
            sys.stderr.write(f'PSIProfile: error parsing results from {search_result}')
            return None

        for query in search_result:
            for hit in list(query):
                for hsp in hit:
                    identity = 1.0 * hsp.ident_num / hsp.aln_span
                    data = [hsp.query.id, hsp.query_start, hsp.query_end, hsp.hit.id, hsp.hit_start, hsp.hit_end,
                            hsp.evalue, identity, str(hsp.aln[0].seq), str(hsp.aln[1].seq)]
                    yield {f: data[i] for i, f in enumerate(PsiProfile.search_result_fields)}

    @staticmethod
    def create_psi_model(seq_id, seq_fasta, profile_database, pssm_build_iterations, pdb_seqres, work_dir, cpus,
                         tmp_dir=None, force=False):
        if not tmp_dir:
            tmp_dir = mkdtemp("structurome_")
        search_result = work_dir + "/profile_search.xml"

        pssm_file = work_dir + "/profile.pssm"

        if force or not os.path.exists(pssm_file):
            PsiProfile.build_profile(seq_fasta, profile_database, pssm_build_iterations, pssm_file, cpus)

        PsiProfile.profile_search(seq_id, pdb_seqres, pssm_file, search_result, cpus)
        search_result = bpsio.parse(search_result, "blast-xml")
        hsps = []
        for query in search_result:
            for hit in list(query):
                for hsp in hit:
                    if hsp.evalue < 10 ** -5:
                        hsps.append(hsp)
        result = Modelome.model_hsps(seq_id, work_dir, hsps, refinement, models_to_generate, assessments, entries, pdb_divided,
                            tmp_dir,max_models=3)

        if not skip_quality:
            for models in result["models"].values():
                for model_path in models:
                    if not os.path.exists(model_path + ".json"):
                        assessment = QMean.assesment(model_path)
                        with open(model_path + ".json", "w") as h:
                            json.dump(assessment, h)





if __name__ == "__main__":
    import argparse
    import sys
    import Bio.SeqIO as bpio
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from SNDG import mkdir
    import fileinput
    from tqdm import tqdm

    def process_pssm(pssm):
        assert os.path.exists(pssm), f'"{pssm}" does not exists'
        if args.alns_dir:
            search_result = f'{args.alns_dir}/{".".join(os.path.basename(pssm).split(".")[:-1])}.xml'
        else:
            search_result = mktemp()
        alns = list(PsiProfile.profile_search(args.database, pssm, search_result, args.cpu))
        if not alns:
            sys.stderr.write(f'{os.path.basename(pssm)} has no hits\n')
        for aln in alns:
            if args.format == "fasta":
                desc = f"{aln['qstart']}-{aln['qend']} identity={aln['identity']} evalue={aln['evalue']}"
                bpio.write(SeqRecord(id=aln["query"], name="", description=desc, seq=Seq(aln['qaln'])), args.output,
                           "fasta")
                bpio.write(SeqRecord(id=aln["hit"], name="", description=f"{aln['hstart']}-{aln['hend']}",
                                     seq=Seq(aln['haln'])), args.output, "fasta")
            else:
                print("\t".join([str(aln[k]) for k in PsiProfile.search_result_fields]))


    parser = argparse.ArgumentParser(description='Profile utils')
    subparsers = parser.add_subparsers(help='commands',  description='valid subcommands',  dest='command')
    pssm = subparsers.add_parser('pssm', help='create PSSM')
    pssm.add_argument('seqs', nargs='?', type=argparse.FileType('r'),
                      default=sys.stdin, help='fasta with AA sequences (default: stdin)')
    pssm.add_argument("-o", '--output', default="./", help='pssm output dir')
    pssm.add_argument("-d", '--database', required=True, help='db used to build the pssm')
    pssm.add_argument("-i", '--iterations', default=3, type=int, help='psi-blast iterations')
    pssm.add_argument('--cpu', default=4, type=int, help='cpus to use')

    search = subparsers.add_parser('search', help='search PSSM')
    search.add_argument('pssms', nargs='*',
                        default="-", help='PSSM files to process')

    search.add_argument("-o", '--output', type=argparse.FileType('w'), default=sys.stdout,
                        help='search result output')
    search.add_argument( '--alns_dir',  default=None,  help='save blast aligments in this folder')
    search.add_argument("-d", '--database', required=True, help='db to be searched with the pssm/s')
    search.add_argument('--cpu', default=4, type=int, help='cpus to use')
    search.add_argument('--format', choices=["table", "fasta"], default="table")

    args = parser.parse_args()

    if args.command == "pssm":
        mkdir(args.output)
        assert os.path.exists(args.output), f"{args.output} could not be created"
        for record in bpio.parse(args.seqs, "fasta"):

            pssm_file = f'{args.output}/{record.id}.pssm'
            query_file = mktemp()
            bpio.write(record,query_file,"fasta")
            if (os.path.exists(pssm_file) and  (os.path.getsize(pssm_file) > 100)):
                print(pssm_file)
            else:
                PsiProfile.build_profile(query_file, args.database, args.iterations, pssm_file, args.cpu)
                if (not os.path.exists(pssm_file)) or  (os.path.getsize(pssm_file) < 100):
                    sys.stderr.write(f'profile for {record.id} could not be generated\n')
                else:
                    print(pssm_file)

        sys.exit(0)

    if args.command == "search":
        if args.format == "table":
            print("\t".join(PsiProfile.search_result_fields))
        if args.alns_dir:
            mkdir(args.alns_dir)
            assert os.path.exists(args.alns_dir), f'{args.alns_dir} could not be created'
        if args.pssms == "-":
            iter_pssm = fileinput.input(args.pssms)
        elif len(args.pssms) > 1:
            iter_pssm = args.pssms
        else:
            with open(args.pssms[0], "r") as file:
                first_line = file.readline()
            if os.path.exists(first_line):
                with open(args.pssms[0], "r") as file:
                    iter_pssm = [x.strip() for x in file.readlines()]
            else:
                iter_pssm = args.pssms

        for pssm in tqdm(iter_pssm) :
            pssm = pssm.strip()
            try:
                process_pssm(pssm)
            except ValueError:
                sys.stderr.write(f'error processing "{pssm}"')

            # ["query", "qstart", "qend", "hit", "hstart", "hend", "evalue", "identity", "qaln", "haln"]

        sys.exit(0)
