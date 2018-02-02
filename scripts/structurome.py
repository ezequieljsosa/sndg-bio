#!/usr/bin/env python

'''
Created on Aug 24, 2015

@author: eze
'''

import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)



import os
import argparse
import multiprocessing
import logging

import Bio.SeqIO as  bpio

from SNDG import init_log
from SNDG.Structure.PsiProfile import create_psi_model

if __name__ == '__main__':



    parser = argparse.ArgumentParser()
    parser.add_argument("-fasta", required=True)
    parser.add_argument("-profile_db", "--profile_db", default="/data/uniprot/uniref/uniref90.fasta")
    parser.add_argument("-pdbs", "--pdb_seqres", default="/data/pdb/processed/seqs_from_pdb.fasta")
    parser.add_argument("-iter", "--pssm_build_iterations", default=3)
    parser.add_argument("-cpus", "--cpus", default=multiprocessing.cpu_count())  # @UndefinedVariable
    parser.add_argument("-o", "--output_dir", default=None)
    parser.add_argument("-l", "--log_path", default=None)
    parser.add_argument("-e", "--entries", default='/data/pdb/entries.idx')
    parser.add_argument("-d", "--pdb_divided", default="/data/pdb/divided/")

    args = parser.parse_args()

    assert os.path.exists(args.profile_db), "%s does not exists" % args.profile_db
    assert os.path.exists(args.pdb_seqres), "%s does not exists" % args.pdb_seqres
    assert os.path.exists(args.entries), "%s does not exists" % args.entries
    assert os.path.exists(args.pdb_divided), "%s does not exists" % args.pdb_divided
    assert os.path.exists(args.fasta), "%s does not exists" % args.pdb_divided

    if not args.output_dir:
        args.output_dir = os.path.dirname(args.fasta)

    if not args.log_path:
        args.log_path = args.fasta + ".log"
    init_log(args.log_path,logging.INFO)
    _log = logging.getLogger("model_fasta")


    # entries_columns =  ["IDCODE", "HEADER", "ACCESSIONDATE", "COMPOUND", "SOURCE", "AUTHORS", "RESOLUTION", "EXPERIMENT"]
    # pdb_entries_df = pd.read_table(args.entries,skiprows=[0,1,2],sep='\t',names=entries_columns)
    def res_fn(y):
        try:
            res = float(y)
        except:
            res = 30
        return res


    entries = {x.split("\t")[0].lower(): res_fn(x.split("\t")[6]) for x in list(open(args.entries))[3:]}

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)



#----------------------



    from collections import defaultdict
    from tqdm import tqdm
    import Bio.SearchIO as bpsio
    from SNDG.Structure.PsiProfile import model_hsps
    con_pdb = [ ]
    good_model = defaultdict(lambda : [])
    def identity(hsp):
        return 1.0 * hsp.ident_num / hsp.aln_span
    for query in bpsio.parse("/data/pext/good_hits.xml","blast-xml"):
        for hit in query:
            if list(hit):
                hsp = list(hit)[0]
                if identity(hsp) >= 0.9:
                    con_pdb.append( hsp.query.id)
                elif identity >= 0.6:
                    good_model[hsp.query.id].append(hsp)


    tuplas = good_model.items()

    with tqdm(tuplas) as pbar:
        for seq,hsps in pbar:
            #pbar.set_description(seq)
            #break
            try:
                model_hsps(seq,"/data/pext_good/",hsps,entries=entries,tmp_dir="/tmp/chain_PDBs")
            except Exception as ex:
                _log.exception(ex)


    import sys
    sys.exit()
#----------------------





    seqs = bpio.parse(args.fasta, "fasta")
    for i, seq in enumerate(seqs):

        try:
            _log.debug((i, seq.id))
            workdir = os.path.abspath(args.output_dir) + "/" + seq.id + "/"
            if not os.path.exists(workdir):
                os.makedirs(workdir)
            fasta = workdir + "seq.fasta"
            bpio.write([seq], fasta, "fasta")

            create_psi_model(seq.id, fasta, args.profile_db, args.pssm_build_iterations, args.pdb_seqres, workdir,
                             args.cpus, entries=entries, pdb_divided=args.pdb_divided)
        except Exception as ex:
            _log.warn(str(ex))
