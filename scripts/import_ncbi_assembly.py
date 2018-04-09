import logging
import multiprocessing
import traceback
from argparse import ArgumentParser, RawDescriptionHelpFormatter

from Bio import Entrez
from peewee import MySQLDatabase

from SNDG import init_log, mkdir
from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
from SNDG.BioMongo.Process.Importer import from_ref_seq, update_proteins, create_proteome
from SNDG.BioMongo.Process.Index import index_seq_collection, build_statistics
from SNDG.BioMongo.Process.Taxon import tax_db
from SNDG.Sequence.ProteinAnnotator import ProteinAnnotator
from SNDG.WebServices.NCBI import ExternalAssembly

Entrez.email = "ezejajaja@hotmail.com"
_log = logging.getLogger(__name__)
if __name__ == "__main__":
    logger = logging.getLogger('peewee')
    logger.setLevel(logging.INFO)
    init_log()

    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-p", "--dbpass", required=True)
    parser.add_argument("-a", "--assemblyAccession", required=True)
    parser.add_argument("-mdb", "--mongodbname", required=True)
    parser.add_argument("-mydbtax", "--mysqldbtaxname", default="bioseqdb")
    parser.add_argument("--cpus", default=multiprocessing.cpu_count())
    parser.add_argument("-mydbunip", "--mysqldbunip", default="unipmap")
    parser.add_argument("-myu", "--mysqldbuser", default="root")

    args = parser.parse_args()
    args.cpus = int(args.cpus)
    mdb = BioMongoDB(args.mongodbname)
    tax_db.initialize(MySQLDatabase(args.mysqldbtaxname, user=args.mysqldbuser, passwd=args.dbpass))
    ProteinAnnotator.connect_to_db(database=args.mysqldbunip, user=args.mysqldbuser, password=args.dbpass)

    assert not mdb.seq_col_exists(args.assemblyAccession), "assembly already exists"
    Entrez.email = "a@a.com"
    assembly_id = Entrez.read(Entrez.esearch(db="assembly", term=args.assemblyAccession , retmax=1))["IdList"][0]
    resource = Entrez.read(Entrez.esummary(db="assembly", id=assembly_id, validate=False))
    try:

        data = resource["DocumentSummarySet"]["DocumentSummary"][0]
        name = data["AssemblyName"]
        genome = str(data["SpeciesName"])
        tax = data["Taxid"]
        status = data["AssemblyStatus"]

        ea = ExternalAssembly(type="assembly", name=name, identifier=assembly_id
                              , assembly_accession=data['AssemblyAccession'], genome=genome
                              , assembly_name=data['AssemblyName'],
                              assembly_status=status, ncbi_tax=int(tax))

        dst_dir = "/data/organismos/" + ea.assembly_accession + "/annotation/"
        mkdir(dst_dir)
        gbpath = ea.download_gbk(dst_dir)
        from_ref_seq(ea.assembly_accession, gbpath, tax=ea.ncbi_tax, tmp_dir=dst_dir)

        tid = int(mdb.db.sequence_collection.find_one({"name": ea.assembly_accession})["tax"]["tid"])
        tmp_dir = "/data/organismos/" + ea.assembly_accession + "/annotation/"
        proteome_dir = "/data/organismos/" + ea.assembly_accession + "/contigs/"
        mkdir(tmp_dir)
        mkdir(proteome_dir)
        protein_fasta = create_proteome(proteome_dir, ea.assembly_accession)
        update_proteins(tmp_dir, protein_fasta, ea.assembly_accession, tid, cpus=args.cpus)
        index_seq_collection(mdb.db, ea.assembly_accession, pathways=False, structure=True)
        build_statistics(mdb.db, ea.assembly_accession)


        jw = JBrowse(db=mdb.db)
        jw.create_genome(ea.assembly_accession)

    except Exception as ex:
        _log.warn(str(ex))
        traceback.print_exc()
        mdb.delete_seq_collection(ea.assembly_accession)
