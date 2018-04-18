#!/usr/bin/env python

'''
Created on Aug 24, 2015

@author: eze
'''
import sys
import tempfile
import warnings

import os
from Bio import BiopythonWarning, BiopythonParserWarning, BiopythonDeprecationWarning, BiopythonExperimentalWarning

warnings.simplefilter('ignore', BiopythonWarning)
warnings.simplefilter('ignore', BiopythonParserWarning)
warnings.simplefilter('ignore', BiopythonDeprecationWarning)
warnings.simplefilter('ignore', BiopythonExperimentalWarning)

import Bio.SeqIO as bpio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
from SNDG.BioMongo.Process.SearchLoader import SearchLoader
from SNDG.Sequence.ProteinAnnotator import ProteinAnnotator
from SNDG.BioMongo.Model.Protein import Protein

import argparse
import logging
from tqdm import tqdm
import Bio.SearchIO as bpsio

from SNDG import init_log, execute

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-o", "--organism", required=True)
    parser.add_argument("-mdb", "--mongo_db", required=True)

    parser.add_argument("-f", "--fasta", default=None)

    parser.add_argument("-mdbp", "--mdb_port",type=int, default=27017)
    parser.add_argument("-mdbh", "--mdb_host", default="127.0.0.1")

    parser.add_argument("-u", "--uniref", default="/data/databases/uniprot/uniref/uniref90.fasta")
    parser.add_argument("-b", "--blast", default=True,
                        help="blastp -qcov_hsp_perc 80 -max_hsps 1 -max_target_seqs 1 -evalue 1e-5 -query {prots} -db {unrefx} -outfmt 6")

    parser.add_argument("-p", "--user_pass", required=True)
    parser.add_argument("-dba", "--db_annotation", default="unipmap")
    parser.add_argument("-dbu", "--user_db", default="root")
    parser.add_argument("-cpu", "--cpu", default="1")

    parser.add_argument("-l", "--log_path", default=None)

    args = parser.parse_args()
    _log = logging.getLogger("protein_annotation")

    if not args.log_path:
        args.log_path = "/tmp/annotation.log"
    init_log(args.log_path, logging.INFO)

    pa = ProteinAnnotator()
    pa.connect_to_db(database=args.db_annotation, user=args.user_db, password=args.user_pass)

    BioMongoDB(args.mongo_db,host=args.mdb_host,port=args.mdb_port)

    if not os.path.exists(args.blast):
        _log.info(args.blast + " does not exists, running blast...")
        if args.fasta:
            assert os.path.exists(args.fasta), args.fasta + " does not exists"
            fasta = args.fasta
        else:
            _log.info("no fasta input, using proteins from the mongo db")
            fasta = tempfile.mktemp()
            with open(fasta, "w") as h:
                for p in Protein.objects(organism=args.organism):
                    r = SeqRecord(id=p.gene[0], description="", seq=Seq(p.seq))
                    bpio.write(r, h, "fasta")


        execute("blastp -qcov_hsp_perc 80 -max_hsps 1 -evalue 1e-5 -query %s -db %s -num_threads %s -outfmt 6 -out %s" %
                (fasta, args.uniref, args.cpu, args.blast))

    with tqdm(list(bpsio.parse(args.blast, "blast-tab"))) as pbar:
        for query in pbar:
            pbar.set_description(query.id)
            try:
                if query[0][0].ident_pct > 0.9:
                    dbxrefs = pa.uniref_annotations(query[0].id, "UniRef90")
                    p = SearchLoader.update_protein_with_dbxref(query.id, dbxrefs, args.organism)
                    p.save()
            except Exception as e:
                _log.warn(e)
    sys.exit(0)
