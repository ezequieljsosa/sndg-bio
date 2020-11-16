import logging
import multiprocessing
import sys

from Bio import Entrez
from tqdm import tqdm

from SNDG import execute, mkdir
from SNDG.WebServices import download_file
from SNDG.WebServices.NCBI import NCBI

Entrez.email = 'A.N.Other@example.com'

_log = logging.getLogger(__name__)

from collections import defaultdict

from SNDG.Annotation.GenebankUtils import GenebankUtils

gut_microbiote_assemblies = [x.strip() for x in """GCA_000712235.1
GCA_002017855.1
GCA_002215605.1
GCF_000144975.1
GCF_000146835.1
GCF_000148995.1
GCF_000151245.1
GCF_000153885.1
GCF_000153905.1
GCF_000153925.1
GCF_000154065.1
GCF_000154085.1
GCF_000154105.1
GCF_000154205.1
GCF_000154285.1
GCF_000154305.1
GCF_000154345.1
GCF_000154365.1
GCF_000154385.1
GCF_000154405.1
GCF_000154425.1
GCF_000154465.1
GCF_000154485.1
GCF_000154505.1
GCF_000154525.1
GCF_000154565.1
GCF_000154805.1
GCF_000154825.1
GCF_000154845.1
GCF_000154865.1
GCF_000154985.1
GCF_000155085.1
GCF_000155205.1
GCF_000155435.1
GCF_000155495.1
GCF_000155835.1
GCF_000155855.1
GCF_000155875.1
GCF_000155955.1
GCF_000155975.1
GCF_000155995.1
GCF_000156015.1
GCF_000156035.2
GCF_000156055.1
GCF_000156075.1
GCF_000156175.1
GCF_000156195.1
GCF_000156215.1
GCF_000156375.1
GCF_000156395.1
GCF_000156495.1
GCF_000156515.1
GCF_000156535.1
GCF_000156655.1
GCF_000157015.1
GCF_000157055.1
GCF_000157115.2
GCF_000157935.1
GCF_000157955.1
GCF_000157975.1
GCF_000157995.1
GCF_000158035.1
GCF_000158055.1
GCF_000158075.1
GCF_000158195.2
GCF_000158315.2
GCF_000158435.2
GCF_000158455.1
GCF_000158475.2
GCF_000158555.2
GCF_000158655.1
GCF_000158835.2
GCF_000159175.1
GCF_000159195.1
GCF_000159215.1
GCF_000159495.1
GCF_000159715.1
GCF_000159915.2
GCF_000159975.2
GCF_000160095.1
GCF_000160175.1
GCF_000160455.2
GCF_000160575.1
GCF_000160595.1
GCF_000161955.2
GCF_000162075.1
GCF_000162115.1
GCF_000162575.1
GCF_000163095.1
GCF_000163735.1
GCF_000163955.1
GCF_000164175.1
GCF_000169015.1
GCF_000169035.1
GCF_000169255.2
GCF_000169475.1
GCF_000172135.1
GCF_000172175.1
GCF_000173355.1
GCF_000173795.1
GCF_000173815.1
GCF_000173975.1
GCF_000174195.1
GCF_000174215.1
GCF_000177015.3
GCF_000178195.1
GCF_000178215.1
GCF_000179075.1
GCF_000185325.1
GCF_000185345.1
GCF_000185665.1
GCF_000185685.2
GCF_000185705.2
GCF_000185845.1
GCF_000186505.1
GCF_000186545.1
GCF_000187265.1
GCF_000187895.1
GCF_000188175.1
GCF_000188195.1
GCF_000191845.1
GCF_000191865.1
GCF_000195635.1
GCF_000204455.1
GCF_000205025.1
GCF_000205165.1
GCF_000213555.1
GCF_000218325.1
GCF_000218405.2
GCF_000220825.1
GCF_000220865.1
GCF_000224635.1
GCF_000224655.1
GCF_000225685.1
GCF_000225705.1
GCF_000225745.1
GCF_000225845.1
GCF_000227195.1
GCF_000227255.2
GCF_000231275.1
GCF_000233455.1
GCF_000233495.1
GCF_000233955.1
GCF_000234155.1
GCF_000234175.1
GCF_000235885.1
GCF_000238035.1
GCF_000238615.1
GCF_000238635.1
GCF_000238655.1
GCF_000238675.1
GCF_000238695.1
GCF_000238735.1
GCF_000238755.1
GCF_000239255.1
GCF_000239295.1
GCF_000239335.1
GCF_000239735.1
GCF_000241405.1
GCF_000242215.1
GCF_000242435.1
GCF_000243175.1
GCF_000243215.1
GCF_000245775.1
GCF_000250875.1
GCF_000261205.1
GCF_000273465.1
GCF_000273585.1
GCF_000296445.1
GCF_000296465.1
GCF_000297815.1
GCF_000315485.1
GCF_000320405.1
GCF_000332875.2
GCF_000345045.1
GCF_000349975.1
GCF_000376405.1
GCF_000381365.1
GCF_000382085.1
GCF_000398925.1
GCF_000411235.1
GCF_000411275.1
GCF_000411295.1
GCF_000411315.1
GCF_000411335.1
GCF_000411415.1
GCF_000412335.1
GCF_000413335.1
GCF_000413355.1
GCF_000413375.1
GCF_000466385.1
GCF_000466465.2
GCF_000466485.1
GCF_000466565.1
GCF_000468015.1
GCF_000469305.1
GCF_000469345.1
GCF_000469445.2
GCF_000479045.1
GCF_000507845.1
GCF_000507865.1
GCF_000517805.1
GCF_000690925.1
GCF_000760655.1
GCF_000763035.1
GCF_000763055.1
GCF_000771165.1
GCF_000969835.1
GCF_000969845.1
GCF_001025135.1
GCF_001025155.1
GCF_001185345.1
GCF_001311295.1
GCF_001315785.1
GCF_001434655.1
GCF_001434945.1
GCF_001435475.1
GCF_001435665.1
GCF_001436305.1
GCF_001941425.1
GCF_002222595.1
GCF_900129655.1
GCF_900167285.1
GCF_001025195.1
GCF_001025215.1
GCF_001434175.1""".split("\n")]

import gzip


class Offtarget(object):
    DEFAULT_GUT_FILENAME = "gut_microbiota.fasta.gz"
    DEFAULT_HUMAN_FILENAME = "human_uniprot100.fa.gz"

    DEG_PROT_URL = {"p": "http://tubic.tju.edu.cn/deg_test/public/download/DEG10.aa.gz",
                    "a": "http://tubic.tju.edu.cn/deg_test/public/download/DEG30.aa.gz",
                    "e": "http://tubic.tju.edu.cn/deg_test/public/download/DEG20.aa.gz"
                    }
    DEG_FAA_NAMES = {
        "a":"degaa-a.dat","p":"degaa-p.dat","e":"degaa-e.dat"
    }

    @staticmethod
    def download_deg(dst="/data/databases/deg/"):
        for x in ["p", "e", "a"]:
            download_file(Offtarget.DEG_PROT_URL[x], f"{dst}/{Offtarget.DEG_FAA_NAMES[x]}.gz", ovewrite=True)
            execute(f"gunzip -f {dst}/{Offtarget.DEG_FAA_NAMES[x]}.gz")
            execute(f"makeblastdb -dbtype prot -in {dst}/{Offtarget.DEG_FAA_NAMES[x]}")

    @staticmethod
    def download_human_prots(dst="/data/databases/human/"):
        file_path = dst + Offtarget.DEFAULT_HUMAN_FILENAME
        unip_url = "https://www.uniprot.org/uniref/?query=uniprot:(taxonomy:%22Homo%20sapiens%20(Human)%20[9606]%22)%20identity:1.0&format=fasta&force=true&compress=yes"
        download_file(unip_url, file_path, ovewrite=True, timeout=120)
        return file_path

    @staticmethod
    def create_human_microbiome(dst="/data/databases/human/", update=False):

        dst_accs = dst + "gut_microbiota_assemblies/"
        mkdir(dst_accs)
        final_file = dst + Offtarget.DEFAULT_GUT_FILENAME

        utils = GenebankUtils()
        with gzip.open(final_file, "wt") as h:
            for accession in tqdm(gut_microbiote_assemblies, file=sys.stderr):
                genome_path = dst_accs + accession + ".genomic.gbff.gz"
                if update or not os.path.exists(genome_path):
                    genome_path = NCBI.download_assembly(accession, dst_accs)
                utils.proteins(genome_path, h)

        return final_file

    @staticmethod
    def count_organism_from_microbiome_blast(tbl_blast_result_path, microbiome_fasta, identity_threshold=0.4,
                                             out_tbl=None):
        prot_org_map = {}
        with (gzip.open(microbiome_fasta, "rt") if microbiome_fasta.endswith(".gz") else open(microbiome_fasta)) as h:
            for line in h:
                if line.startswith(">"):
                    seqid = line.split()[0].strip().replace(">", "")
                    try:
                        org = line.replace("[[", "[").split("[")[1].strip()[:-1]
                    except IndexError:
                        err = "fasta does not have the organism name at the fasta header."
                        err += "example: >HMPREF1002_RS00015 alpha/beta hydrolase [Porphyromonas sp. 31_2]"
                        raise LookupError(err)

                    prot_org_map[seqid] = org

        query_orgs = defaultdict(lambda: [])
        with open(tbl_blast_result_path) as h:
            for l in list(h)[1:]:
                query, hit, identity = l.split()[:3]
                identity = float(identity) / 100.0
                if identity_threshold <= identity:
                    query_orgs[query].append(prot_org_map[hit])
        for query, hits in query_orgs.items():
            query_orgs[query] = set(hits)
        if out_tbl:
            with open(out_tbl, "w") as h:
                for query, hits in query_orgs.items():
                    h.write("\t".join([query, str(len(hits)), ";".join(hits)]) + "\n")
        return query_orgs

    @staticmethod
    def offtargets(proteome, dst_resutls, offtarget_db, cpus=multiprocessing.cpu_count()):
        cmd = f"blastp -evalue 1e-5 -max_hsps 1 -outfmt 6  -db {offtarget_db} -query {proteome} -out {dst_resutls} -num_threads {cpus}|awk '$3>50'"
        execute(cmd)
        return dst_resutls


if __name__ == "__main__":
    from SNDG import init_log
    import argparse
    import os
    from SNDG.Sequence import smart_parse

    parser = argparse.ArgumentParser(description='Offtarget Utilities')

    subparsers = parser.add_subparsers(help='commands', description='valid subcommands', required=True, dest='command')

    gut_download = subparsers.add_parser('download', help='Download offtarget data')
    gut_download.add_argument('-db', '--databases', choices=["all","deg", "human", "gut_microbiote"], default="all")
    gut_download.add_argument('-o', '--output', help="output_directory", default="/data/databases/")
    gut_download.add_argument('--force', action="store_true")

    gut_microbiote_blast = subparsers.add_parser('gut_microbiote_blast',
                                                 help='Runs blastp against gut microbiote and counts organisms')
    gut_microbiote_blast.add_argument('input_faa')
    gut_microbiote_blast.add_argument('-o', '--output', help="output_directory", default="./")
    gut_microbiote_blast.add_argument('-db', '--database', help="gut microbiome fasta",
                                      default="/data/databases/human/gut_microbiota.fasta.gz")
    gut_microbiote_blast.add_argument('--cpus', default=multiprocessing.cpu_count())
    gut_microbiote_blast.add_argument('--force', action="store_true")

    args = parser.parse_args()

    init_log()

    if args.command == "download":
        if args.databases in ["all", "gut_microbiote"]:
            path = f'{args.output}/gut_microbiote/{Offtarget.DEFAULT_GUT_FILENAME}'
            if args.force or not os.path.exists(path):
                path = Offtarget.create_human_microbiome(dst=args.output)
            else:
                sys.stderr.write(f'{path} already exists, overwrite using --force')

            filename = os.path.basename(path)
            execute(f"zcat {path} | makeblastdb -title {filename} -out  {args.output}/{filename} -dbtype prot -in -")
        if args.databases in ["all", "human"]:
            path = f'{args.output}/human/{Offtarget.DEFAULT_HUMAN_FILENAME}'
            if args.force or not os.path.exists(path):
                path = Offtarget.download_human_prots(dst=args.output)
            else:
                sys.stderr.write(f'{path} already exists, overwrite using --force')

            filename = os.path.basename(path)
            execute(f"zcat {path} | makeblastdb -title {filename} -out  {args.output}/{filename} -dbtype prot -in -")
        if args.databases in ["all", "deg"]:
            mkdir(f'{args.output}/deg/')
            Offtarget.download_deg(f'{args.output}/deg/')
    elif args.command == "gut_microbiote_blast":
        blast_gut_path = f'{args.output}/gut_microbiome.blast.tbl'
        gut_result_path = f'{args.output}/gut_microbiome.tbl'
        if not os.path.exists(args.database + ".phr"):
            raise FileNotFoundError(f"{args.database} index files could not be found. Run makeblastdb")
        if args.force or not os.path.exists(blast_gut_path):
            Offtarget.offtargets(args.input_faa, blast_gut_path, offtarget_db=args.database, cpus=args.cpus)
        else:
            sys.stderr.write(f'{blast_gut_path} already exists, overwrite using --force')

        Offtarget.count_organism_from_microbiome_blast(blast_gut_path, args.database, identity_threshold=0.5,
                                                       out_tbl=gut_result_path)
