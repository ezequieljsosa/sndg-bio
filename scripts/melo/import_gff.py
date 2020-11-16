import logging
from SNDG import Struct, init_log, mkdir

init_log()
import glob
import pymongo
from SNDG.WebServices.NCBI import NCBI
from SNDG.WebServices.Offtargeting import Offtargeting
import Bio.SeqIO as bpio
from tqdm import tqdm

logging.getLogger("peewee").setLevel(logging.WARN)
from peewee import MySQLDatabase
from SNDG.BioMongo.Process.Taxon import tax_db
from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
from SNDG.Sequence.ProteinAnnotator import ProteinAnnotator, Mapping
from SNDG.BioMongo.Process.Importer import from_ref_seq, update_proteins, import_prop_blast
from SNDG.BioMongo.Process.BioDocFactory import BioDocFactory
from SNDG.BioMongo.Model.Protein import Protein
from SNDG.Network.KEGG import Kegg
from SNDG.BioMongo.Process.Importer import _common_annotations, _protein_iter, import_kegg_annotation, \
    index_seq_collection, build_statistics, load_pathways
from BCBio import GFF
from SNDG.BioMongo.Process.Taxon import Tax

tax_db.initialize(MySQLDatabase('bioseqdb', user='root', passwd="mito"))
mdb = BioMongoDB("tdr")
mysqldb = ProteinAnnotator.connect_to_db(database="unipmap", user="root", password="mito")

name = "MeloI"
org = "Meloidogyne incognita"
gff = "/mnt/Data/data/projects/melo/target/Meloidogyne_incognita_V3_20131230.gff3"
fna = "/mnt/Data/data/projects/melo/target/Meloidogyne_incognita_V3_20131230.fna"
protsfasta = "/mnt/Data/data/projects/melo/target/meloidogyne_incognita.v3.2015.fpa"

#
# seqs = {r.id: r.seq for r in bpio.parse(fna, "fasta")}
# protseqs = {r.id: str(r.seq) for r in bpio.parse(protsfasta, "fasta")}
#
#
# seqCol = BioDocFactory.create_genome(name, Struct(description=org,annotations={"organism":org}), 6306, Tax)
# seqCol.description = org + " (https://meloidogyne.inra.fr/Downloads/Meloidogyne-incognita-V2-2017)"
# seqCol.organism = org
# seqCol.save()
#
# tmp_dir = "/tmp/" + name + "/"
# mkdir(tmp_dir)
# gene_ids = {}
#
# ann = list(GFF.parse(gff))
# with tqdm(ann) as pbar:
#     for contig in pbar:
#         pbar.set_description(contig.id)
#         if len(contig.seq) > 15000000:
#             contig.seq = ""
#         contigDoc, gene_ids2 = BioDocFactory.create_contig(
#             contig, seqCol, type_map=
#             {"rRNA": "rRNA", "ncRNA": "ncRNA", "mRNA": "gene", "gene": "gene",
#              NCBI.f_CDS: NCBI.f_CDS, "rRNA": "rRNA", "tRNA": "tRNA", "tmRNA": "tmRNA"},
#             extract_annotation_feature=lambda x:x.sub_features[0],
#         )
#         gene_ids.update(gene_ids2)
#         contigDoc.save()
#
# prots = []
# with tqdm(ann) as pbar:
#     for contig in pbar:
#         for f in contig.features:
#             if len(f.sub_features) and f.sub_features[0].type == "mRNA":
#                 locus_tag = f.qualifiers["Name"][0]
#                 cds_f = f.sub_features[0]
#                 seq = protseqs[ locus_tag ]
#                 protDoc = Protein(seq=seq, name=locus_tag,gene=[locus_tag],alias=[locus_tag],description="",
#                             ontologies=[],properties=[])
#                 protDoc.gene_id = gene_ids["mRNA:" + locus_tag]
#                 protDoc.organism = name
#                 protDoc.auth = str(BioMongoDB.demo_id)
#                 protDoc.seq_collection_id = seqCol
#                 prots.append(protDoc)
#                 if pbar.n and ((pbar.n % 1000) == 0):
#                     Protein.objects.insert(prots)
#                     prots = []
#         if prots:
#             Protein.objects.insert(prots)
#
# _common_annotations(name, tmp_dir, cpu=1)


# mdb.protein_fasta("/tmp/" + name + "/genome.fasta", name)
# update_proteins("/tmp/" + name + "/", "/tmp/" + name + "/genome.fasta", name, 6306, db_init=mysqldb)


# Offtargeting.offtargets("/tmp/" + name + "/genome.fasta",
#                         "/data/organismos/" + name + "/annotation/offtarget/",
#                         offtarget_dbs=[  "/data/databases/deg/degaa-e.dat",
#                                          "/data/databases/human/gencode.v17.pc_translations.fa",
#                                          "/data/databases/human/gut_microbiota.fasta"]
#                         )

# import_prop_blast(mdb.db, name, "hit_in_deg",
#                   "/data/organismos/" + name + "/annotation/offtarget/degaa-e.tbl",
#                   "table", "Hit in DEG database",
#                   value_fn=lambda x: x.identity > 70,
#                   default_value=True,
#                   no_hit_value=False, choices=["True", "False"], type="value", defaultOperation="equal")
#
# import_prop_blast(mdb.db, name, "human_offtarget",
#                   "/data/organismos/" + name + "/annotation/offtarget/gencode.tbl",
#                   "table", "Human offtarget score (1 - best hit identity)",
#                   value_fn=lambda x: 1 - (x.identity * 1.0 / 100),
#                   default_value=0.4,
#                   no_hit_value=1)
# import_prop_blast(mdb.db, name, "gut_microbiota_offtarget",
#                   "/data/organismos/" + name + "/annotation/offtarget/gut_microbiota.tbl",
#                   "table", "Gut microbiota offtarget score (1 - best hit identity)",
#                   value_fn=lambda x: 1 - (x.identity * 1.0 / 100),
#                   default_value=0.4,
#                   no_hit_value=1)

# kegg = Kegg(ko_list_path="/data/databases/kegg/ko.txt",
#             brittle_pw_path="/data/databases/kegg/ko/",
#             kgmls_dir="/data/databases/kegg/ko/")
# kegg.init()
# ilex_data = "/data/organismos/MeloI/annotation/query.ko"
# kegg.read_annotation(ilex_data)
#
# import_kegg_annotation(mdb.db,"MeloI",kegg)

load_pathways(name, "/data/organismos/MeloI/annotation/pwtools/allreactions.sbml", mdb.db,
              "/data/organismos/MeloI/annotation/pwtools/",
              gregexp="\(([\-\w\.]+)\)", filter_file="allfilters_con_c.dat")

# index_seq_collection(mdb.db,name,pathways=False,go=True,keywords=True,ec=True,organism_idx=True,structure=True)
# build_statistics(mdb.db,name)
