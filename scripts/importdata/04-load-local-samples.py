import logging
import traceback

from Bio import Entrez
from tqdm import tqdm

from SNDG import init_log, mkdir
from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
from SNDG.BioMongo.Process.Importer import from_ref_seq, update_proteins,create_proteome
from SNDG.BioMongo.Process.Taxon import tax_db
from SNDG.WebServices.NCBI import ExternalAssembly,mysql_db
from peewee import MySQLDatabase
from SNDG.Sequence.ProteinAnnotator import ProteinAnnotator, Mapping
from SNDG.BioMongo.Process.Index import index_seq_collection, build_statistics

Entrez.email = "ezejajaja@hotmail.com"
_log = logging.getLogger(__name__)
if __name__ == "__main__":
    logger = logging.getLogger('peewee')
    logger.setLevel(logging.INFO)
    init_log()


    mdb = BioMongoDB("saureus")
    tax_db.initialize(MySQLDatabase('bioseqdb', user='root', passwd="mito"))
    mysql_db.initialize(MySQLDatabase('sndg', user='root', passwd="mito"))
    assemblies = list(ExternalAssembly.select().where(ExternalAssembly.sample_source.is_null(False)))

    ProteinAnnotator.connect_to_db(database="unipmap", user="root", password="mito")
    with tqdm(assemblies) as pbar:
        for x in pbar:
            if mdb.seq_col_exists(x.assembly_accession):
                continue
            pbar.set_description(x.assembly_accession)
            try:
                dst_dir = "/data/organismos/" + x.assembly_accession + "/annotation/"
                mkdir(dst_dir)
                gbpath = x.download_gbk(dst_dir)
                from_ref_seq(x.assembly_accession, gbpath, tax=x.ncbi_tax, tmp_dir=dst_dir)

                tid = int(mdb.db.sequence_collection.find_one({"name": x.assembly_accession})["tax"]["tid"])
                tmp_dir = "/data/organismos/" + x.assembly_accession + "/annotation/"
                proteome_dir = "/data/organismos/" + x.assembly_accession + "/contigs/"
                mkdir(tmp_dir)
                mkdir(proteome_dir)
                protein_fasta = create_proteome(proteome_dir, x.assembly_accession)
                update_proteins(tmp_dir, protein_fasta, x.assembly_accession, tid,cpus=3)
                index_seq_collection(mdb.db, x.assembly_accession, pathways=False, structure=True)
                build_statistics(mdb.db, x.assembly_accession)
            except Exception as ex:
                _log.warn(str(ex))
                traceback.print_exc()
                mdb.delete_seq_collection(x.assembly_accession)

