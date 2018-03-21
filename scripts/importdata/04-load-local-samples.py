import logging

import xmltodict
from Bio import Entrez
from SNDG import init_log, mkdir
from SNDG.WebServices.NCBI import ExternalAssembly
from tqdm import tqdm
from SNDG.BioMongo.Process.Importer import from_ref_seq,update_proteins
from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
import traceback

Entrez.email = "ezejajaja@hotmail.com"

if __name__ == "__main__":
    logger = logging.getLogger('peewee')
    logger.setLevel(logging.INFO)
    init_log()
    assemblies = list(ExternalAssembly.select().where(ExternalAssembly.sample_source.is_null(False)))

    mdb = BioMongoDB("saureus")

    for x in tqdm(assemblies):
        try:
            dst_dir = "/data/organismos/" + x.assembly_accession + "/annotation/"
            mkdir(dst_dir)
            gbpath = x.download_gbk(dst_dir)
            from_ref_seq(x.assembly_accession, gbpath, tax=x.ncbi_tax, tmp_dir=dst_dir)

            tid = int(mdb.db.sequence_collection.find_one({"name": seq_col_name})["tax"]["tid"])
            tmp_dir = "/data/organismos/" + x.assembly_accession + "/annotation/"
            proteome_dir = "/data/organismos/" + x.assembly_accession + "/contigs/"
            mkdir(tmp_dir)
            mkdir(proteome_dir)
            protein_fasta = create_proteome(proteome_dir, seq_col_name)
            update_proteins(tmp_dir, protein_fasta, seq_col_name, tid)
            index_seq_collection(mdb.db, x.assembly_accession, pathways=False, structure=True)
            build_statistics(mdb.db, x.assembly_accession)
        except Exception as ex:
            _log.warn(str(ex))
            traceback.print_exception(*exc_info)
