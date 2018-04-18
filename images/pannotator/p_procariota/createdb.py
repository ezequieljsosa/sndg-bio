#!/usr/bin/env python

import os
import time
import shutil

from SNDG import mkdir, execute, execute_from, init_log
from SNDG.WebServices import download_file

init_log("/tmp/createdb.log")


def old_or_inexistent(filepath, period=30):
    return not os.path.exists(filepath) or (((time.time() - os.path.getatime(filepath)) / 60 / 60 / 24) > period)


os.environ["http_proxy"] = "http://proxy.fcen.uba.ar:8080"
os.environ["ftp_proxy"] = "http://proxy.fcen.uba.ar:8080"



if not os.path.exists("/data/cog/whog"):
    mkdir("/data/cog/")
    download_file("ftp://ftp.ncbi.nih.gov/pub/COG/COG/whog",
                  "/data/cog/whog")

if not os.path.exists("/data/cog/myva"):
    mkdir("/data/cog/")
    download_file("ftp://ftp.ncbi.nih.gov/pub/COG/COG/myva",
                  "/data/cog/myva")
    execute("formatdb -i /data/cog/myva -o T")

if not os.path.exists("/data/ec/PRIAM_MAR15/priam"):
    mkdir("/data/ec/")
    download_file("http://priam.prabi.fr/REL_MAR15/Distribution.zip",
                  "/data/ec/PRIAM_MAR15.zip")
    execute_from("unzip /data/ec/PRIAM_MAR15.zip; exit 0;", "/data/ec/",retcodes=[0,1])

    execute_from("ls /data/ec/PRIAM_MAR15/PROFILES/*.chk > priam", "/data/ec/PRIAM_MAR15/")
    execute_from("formatrpsdb -i /data/ec/PRIAM_MAR15/priam -o T", "/data/ec/PRIAM_MAR15/")

if not os.path.exists("/data/pfamtigrfam/tirgfam.hmm"):
    mkdir("/data/pfamtigrfam/INFO")
    download_file("ftp://ftp.jcvi.org/pub/data/TIGRFAMs/TIGRFAMs_15.0_HMM.LIB.gz",
                  "/data/pfamtigrfam/TIGRFAMs_15.0_HMM.LIB.gz")
    execute("gunzip /data/pfamtigrfam/TIGRFAMs_15.0_HMM.LIB.gz")
    execute_from("hmmconvert TIGRFAMs_15.0_HMM.LIB > tirgfam.hmm", "/data/pfamtigrfam/")

    download_file("ftp://ftp.jcvi.org/pub/data/TIGRFAMs/TIGRFAMs_15.0_INFO.tar.gz",
                  "/data/pfamtigrfam/INFO/TIGRFAMs_15.0_INFO.tar.gz")
    execute_from("tar -xvf /data/pfamtigrfam/INFO/TIGRFAMs_15.0_INFO.tar.gz","/data/pfamtigrfam/INFO/")
    os.remove("/data/pfamtigrfam/INFO/TIGRFAMs_15.0_INFO.tar.gz")

if not os.path.exists("/data/pfamtigrfam/pfam2go.txt"):
    download_file("http://geneontology.org/external2go/pfam2go", "/data/pfamtigrfam/pfam2go.txt")
if not os.path.exists("/data/pfamtigrfam/TIGRFAMS_GO_LINK"):
    download_file("ftp://ftp.jcvi.org/pub/data/TIGRFAMs/TIGRFAMS_GO_LINK", "/data/pfamtigrfam/TIGRFAMS_GO_LINK")

if old_or_inexistent("/data/pfamtigrfam/Pfam-A.hmm", 150):
    download_file("ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz",
              "/data/pfamtigrfam/Pfam-A.hmm.gz")
    execute("gunzip /data/pfamtigrfam/Pfam-A.hmm.gz")

if old_or_inexistent("/data/pfamtigrfam/pfamatigrfam.hmm", 150):
    execute_from("cat tirgfam.hmm Pfam-A.hmm > pfamatigrfam.hmm", "/data/pfamtigrfam/")
    execute("hmmpress /data/pfamtigrfam/pfamatigrfam.hmm")
    shutil.copy("/app/p_procariota/shorthmm_generator.sh", "/data/pfamtigrfam/")
    execute_from("./shorthmm_generator.sh", "/data/pfamtigrfam/")

if old_or_inexistent("/data/uniprot/goa/"):
    mkdir("/data/uniprot/goa/")
    download_file("ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/goa_uniprot_all.gpa.gz",
                  "/data/uniprot/goa/goa_uniprot_all.gpa.gz", ovewrite=True)
    shutil.copy("/app/p_procariota/goagenerator.py", "/data/uniprot/goa/")
    execute("gunzip /data/uniprot/goa/goa_uniprot_all.gpa.gz")
    execute_from("./goagenerator.py", "/data/uniprot/goa/",retcodes=[0,1])

if old_or_inexistent("/data/uniprot/gsm/uniprotheaders",150):
    mkdir("/data/uniprot/gsm")

    if old_or_inexistent("/data/uniprot/gsm/uniprot_sprot.fasta.gz"):
        download_file(
            "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
            "/data/uniprot/gsm/uniprot_sprot.fasta.gz", ovewrite=True)
        try:
            execute("gunzip /data/uniprot/gsm/uniprot_sprot.fasta.gz")
        except:
            os.remove("/data/uniprot/gsm/uniprot_sprot.fasta.gz")

    if old_or_inexistent("/data/uniprot/gsm/uniprot_trembl.fasta.gz"):
        download_file(
            "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz",
        "/data/uniprot/gsm/uniprot_trembl.fasta.gz", ovewrite=True)
        try:
            execute("gunzip /data/uniprot/gsm/uniprot_trembl.fasta.gz")
        except:
            os.remove("/data/uniprot/gsm/uniprot_trembl.fasta.gz")

    shutil.copy("/app/p_procariota/genesymbolmapgenerator.sh", "/data/uniprot/gsm/")
    shutil.copy("/app/p_procariota/genesymbolmapgenerator.py", "/data/uniprot/gsm/")
    execute_from("./genesymbolmapgenerator.sh", "/data/uniprot/gsm/",retcodes=[0,1])

if old_or_inexistent("/data/uniprot/uniref/uniref50/uniref50.fasta"):
    mkdir("/data/uniprot/uniref/uniref50")
    download_file("ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz",
                  "/data/uniprot/uniref/uniref50/uniref50.fasta.gz", ovewrite=True)
    execute("gunzip /data/uniprot/uniref/uniref50/uniref50.fasta.gz")

if old_or_inexistent("/data/uniprot/uniref/uniref50/uniref50.fasta.pal"):
    execute("makeblastdb -dbtype prot -in /data/uniprot/uniref/uniref50/uniref50.fasta")
