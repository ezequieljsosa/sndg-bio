#!/usr/bin/env python

'''
Created on Aug 24, 2015

@author: eze
'''

import warnings

from Bio import BiopythonWarning, BiopythonParserWarning

warnings.simplefilter('ignore', BiopythonWarning)
warnings.simplefilter('ignore', BiopythonParserWarning)

import os
import argparse
import multiprocessing
import logging
from tqdm import tqdm

import Bio.SeqIO as  bpio

from SNDG import init_log
from SNDG.Structure.PsiProfile import PsiProfile
from SNDG.Structure.PDBs import PDBs

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-fasta", required=True)
    parser.add_argument("-profile_db", "--profile_db", default="/data/uniprot/uniref/uniref90/uniref90.fasta")
    parser.add_argument("-pdbs", "--pdb_seqres", default="/data/pdb/processed/seqs_from_pdb.fasta")
    parser.add_argument("-iter", "--pssm_build_iterations", default=3)
    parser.add_argument("-cpus", "--cpus", default=multiprocessing.cpu_count())  # @UndefinedVariable
    parser.add_argument("-o", "--output_dir", default=None)
    parser.add_argument("-l", "--log_path", default=None)
    parser.add_argument("-e", "--entries", default='/data/pdb/entries.idx')
    parser.add_argument("-d", "--pdb_divided", default="/data/pdb/divided/")
    parser.add_argument("--skipProfile", action='store_true', default=False)
    parser.add_argument("--skipQuality", action='store_true', default=False)


    args = parser.parse_args()

    assert ((not args.skipProfile) or os.path.exists(args.profile_db + ".00.phr"),
            "%s does not exists or it is not indexed" % args.profile_db)
    assert os.path.exists(args.pdb_seqres), "%s does not exists" % args.pdb_seqres
    assert os.path.exists(args.entries), "%s does not exists" % args.entries
    assert os.path.exists(args.pdb_divided), "%s does not exists" % args.pdb_divided
    assert os.path.exists(args.fasta), "%s does not exists" % args.fasta

    if not args.output_dir:
        args.output_dir = os.path.dirname(args.fasta)

    if not args.log_path:
        args.log_path = args.fasta + ".log"
    init_log(args.log_path, logging.INFO)
    _log = logging.getLogger("model_fasta")


    def res_fn(y):
        try:
            res = float(y)
        except:
            res = 30
        return res


    entries = {x.split("\t")[0].lower(): res_fn(x.split("\t")[6]) for x in list(open(args.entries))[3:]}

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    seqs = list(bpio.parse(args.fasta, "fasta"))
    with tqdm(seqs) as pbar:
        for seq in pbar:
            try:
                pbar.set_description(seq.id)
                workdir = os.path.abspath(args.output_dir) + "/" + seq.id + "/"
                if not os.path.exists(workdir):
                    os.makedirs(workdir)
                fasta = workdir + "seq.fasta"
                bpio.write([seq], fasta, "fasta")

                PsiProfile.create_psi_model(seq.id, fasta, args.profile_db,
                                            args.pssm_build_iterations, args.pdb_seqres, workdir,
                                            args.cpus, entries=entries, pdb_divided=args.pdb_divided,
                                            skip_profile=args.skipProfile,skip_quality=args.skipQuality)
            except Exception as ex:
                _log.warning(str(ex))
