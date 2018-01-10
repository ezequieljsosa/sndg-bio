#!/usr/bin/env python

import os
import time

from SNDG import mkdir, execute, execute_from, init_log
from SNDG.WebServices import download_file
from SNDG.Structure.PDBs import PDBs

init_log("/tmp/createdb.log")

def old_or_inexistent(filepath, period=30):
    return not os.path.exists(filepath) or (((time.time() - os.path.getatime(filepath)) / 60 / 60 / 24) > period)


os.environ["http_proxy"] = "http://proxy.fcen.uba.ar:8080"
os.environ["ftp_proxy"] = "http://proxy.fcen.uba.ar:8080"

mkdir("/data/pdb/")
download_file("ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/entries.idx", "/data/pdb/entries.idx")

pdbs = PDBs("/data/pdb/")
pdbs.download_pdb_seq_ses()
pdbs.update_pdb_dir()
pdbs.pdbs_seq_for_modelling()
execute("makeblastdb -dbtype prot -in /data/pdb/processed/seqs_from_pdb.fasta")

if old_or_inexistent("/data/uniprot/uniref/uniref90.fasta"):
    mkdir("/data/uniprot/uniref")
    download_file("ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz",
                  "/data/uniprot/uniref90.fasta.gz", ovewrite=True)
    execute("gunzip /data/uniprot/uniref90.fasta.gz")

if old_or_inexistent("/data/uniprot/uniref/uniref90.fasta.pal"):
    execute("makeblastdb -dbtype prot -in /data/uniprot/uniref90.fasta")
