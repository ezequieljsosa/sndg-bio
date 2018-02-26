"""

"""
import os

import sys
from tqdm import tqdm
import subprocess
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.SeqIO as bpio
from SNDG import init_log
from SNDG.BioMongo.Model.Protein import Protein
from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB

from SNDG import mkdir
from SNDG.BioMongo.Process.SearchLoader import load_hmm, load_blast_pdb
from SNDG.Sequence import smart_parse as sp
from SNDG.BioMongo.Process.BioDocFactory import BioDocFactory
from SNDG.BioMongo.Process.Taxon import Tax
from SNDG.Sequence.Hmmer import Hmmer
from SNDG.Sequence.so import SO_TERMS
from SNDG.WebServices.NCBI import NCBI


def _protein_iter(contigIterator):
    for c in contigIterator:
        for f in c.features:
            if f.type == "CDS":
                if "translation" in f.qualifiers:
                    yield (SeqRecord(id=f.qualifiers["locus_tag"][0],
                                     description="",
                                     annotations=f.qualifiers,
                                     seq=Seq(f.qualifiers["translation"][0])), f,)


def tritryp_protein_iter(contigIterator):
    for c in contigIterator:
        for f in c.features:
            if f.type == "gene":
                for fs in f.sub_features:
                    if fs.type == "mRNA":
                        cds = [x for x in fs.sub_features if x.type == "CDS"][0]
                        seq = cds.extract(c).translate().seq
                        yield (SeqRecord(id=f.qualifiers["ID"][0],
                                         description="",
                                         annotations=f.qualifiers,
                                         seq=seq), fs,)


def from_ref_seq(name, gbpath, tax=None, tmp_dir=None):
    for contig in sp(gbpath):
        seqCol = BioDocFactory.create_genome(name, contig, tax, Tax)
        seqCol.save()
        break
    if not tmp_dir:
        tmp_dir = "/tmp/" + name + "/"
    mkdir(tmp_dir)
    gene_ids = {}
    with tqdm(sp(gbpath)) as pbar:
        for contig in pbar:
            pbar.set_description(contig.id)
            if len(contig.seq) > 15000000:
                contig.seq = ""
            contigDoc, gene_ids2 = BioDocFactory.create_contig(contig, seqCol, type_map=
            {"rRNA": "rRNA", "ncRNA": "ncRNA", NCBI.f_mRNA: NCBI.f_mRNA, "gene": "gene",
             NCBI.f_CDS: NCBI.f_CDS, "rRNA": "rRNA", "tRNA": "tRNA", "tmRNA": "tmRNA"})
            gene_ids.update(gene_ids2)
            contigDoc.save()

    prots = []

    with tqdm(_protein_iter(sp(gbpath))) as pbar:
        for (protein, cds_f) in pbar:
            if "locus_tag" in cds_f.qualifiers:
                protDoc = BioDocFactory.create_protein(protein, cds_f)
                if len(protDoc.seq) > 30000:
                    raise Exception("No existen proteinas tan largas...")
                protDoc.gene_id = gene_ids[cds_f.qualifiers["locus_tag"][0]]
                protDoc.organism = name
                protDoc.auth = BioMongoDB.demo_id
                protDoc.seq_collection_id = seqCol
                prots.append(protDoc)
                if pbar.n and ((pbar.n % 1000) == 0):
                    Protein.objects.insert(prots)
                    prots = []
    if prots:
        Protein.objects.insert(prots)

    _common_annotations(name,tmp_dir)


def from_TriTrypDB(name, gff, fasta, tax, tmp_dir=None):
    # genome = {x.id: x for x in sp(fasta)}
    # from BCBio import GFF
    # import re
    # annotation =list( GFF.parse(gff, base_dict=genome))
    # contig = annotation[0]
    #
    # seqCol = BioDocFactory.create_genome(name, contig, tax, Tax)
    # seqCol.save()
    #
    # if not tmp_dir:
    #     tmp_dir = "/tmp/" + name + "/"
    # mkdir(tmp_dir)
    # gene_ids = {}
    # with tqdm(annotation) as pbar:
    #     for contig in pbar:
    #         pbar.set_description(contig.id)
    #         if len(contig.seq) > 15000000:
    #             contig.seq = ""
    #         contigDoc, gene_ids2 = BioDocFactory.create_contig(contig, seqCol, type_map=
    #         {"rRNA": "rRNA", "ncRNA": "ncRNA", NCBI.f_mRNA: "gene", "exon": "exon", "gene": "gene",
    #          NCBI.f_CDS: NCBI.f_CDS, "rRNA": "rRNA", "tRNA": "tRNA", "tmRNA": "tmRNA",
    #          "snoRNA": "snoRNA",
    #          "three_prime_UTR": "three_prime_UTR", "five_prime_UTR": "five_prime_UTR"
    #          })
    #         gene_ids.update(gene_ids2)
    #         contigDoc.save()
    # prots = []
    # with tqdm(tritryp_protein_iter(annotation)) as pbar:
    #     for (protein, cds_f) in pbar:
    #
    #         protDoc = Protein(seq=str(protein.seq), name=protein.id)
    #
    #         if "description" in cds_f.qualifiers:
    #             protein_description = cds_f.qualifiers['description'][0]
    #         elif "Note" in cds_f.qualifiers:
    #             protein_description = cds_f.qualifiers['Note'][0]
    #         elif "product" in cds_f.qualifiers:
    #             protein_description = cds_f.qualifiers['product'][0]
    #         else:
    #             protein_description = ""
    #
    #         protDoc.description = protein_description
    #
    #         gos = []
    #         if "Ontology_term" in cds_f.qualifiers:
    #             gos = [x.lower() for x in cds_f.qualifiers["Ontology_term"] if
    #                    "GO:" in x and (x not in ["GO:0008150", "GO:0003674", "GO:0005575"])]
    #
    #         note = cds_f.qualifiers["Note"][0].split(" ")[0] if "Note" in cds_f.qualifiers else ""
    #         ecs = ["ec:" + note] if re.match('^[0-9]+\.[0-9\-]+\.[0-9\-]+\.[0-9\-]$', note) else []
    #         ontologies = list(set(ecs + gos))
    #
    #         protDoc.gene = [protein.id]
    #         protDoc.ontologies = ontologies
    #         protDoc.alias = [protein.id]
    #
    #         if len(protDoc.seq) > 30000:
    #             raise Exception("No existen proteinas tan largas...")
    #         protDoc.gene_id = gene_ids[protein.id]
    #         protDoc.organism = name
    #         protDoc.auth = BioMongoDB.demo_id
    #         protDoc.seq_collection_id = seqCol
    #         prots.append(protDoc)
    #         if pbar.n and ((pbar.n % 1000) == 0):
    #             Protein.objects.insert(prots)
    #             prots = []
    # if prots:
    #     Protein.objects.insert(prots)
    #
    _common_annotations(name,tmp_dir)

def _common_annotations(name,tmp_dir):
    protein_fasta = tmp_dir + "/proteins.fasta"
    if not os.path.exists(protein_fasta) or (not os.path.getsize(protein_fasta)):
         with open(protein_fasta, "w") as h:
             for p in Protein.objects(organism=name).no_cache():
                 bpio.write(SeqRecord(id=p.gene[0], description="", seq=Seq(p.seq)), h, "fasta")

    if not Protein.objects(
            __raw__={"organism": name, "features.type": SO_TERMS["polypeptide_structural_motif"]}).count():
        blast_result = tmp_dir + "/pdb_blast.xml"
        pdbs_path = "/data/databases/pdb/processed/seqs_from_pdb.fasta"
        if not os.path.exists(blast_result):
            cmd = "blastp -qcov_hsp_perc 80 -max_hsps 1 -evalue 1e-5 -query %s -db %s -num_threads 1 -outfmt 5 -out %s"
            subprocess.call(cmd % (protein_fasta, pdbs_path, blast_result), shell=True)

        load_blast_pdb(name, blast_result)
        #os.remove(blast_result)

    if not Protein.objects(__raw__={"organism": name, "features.type": SO_TERMS["polypeptide_domain"]}).count():
        hmm_result = tmp_dir + "/domains.hmm"
        Hmmer(protein_fasta, output_file=hmm_result).query()
        load_hmm(name, hmm_result)
        os.remove(hmm_result)


if __name__ == '__main__':
    init_log()

    import logging

    logging.getLogger("peewee").setLevel(logging.WARN)
    from peewee import MySQLDatabase
    from SNDG.BioMongo.Process.Taxon import tax_db

    tax_db.initialize(MySQLDatabase('bioseqdb', user='root', passwd="mito"))
    mdb = BioMongoDB("tdr")
    # mdb.delete_seq_collection("GCF_001624625.1")
    # from_ref_seq("GCF_001624625.1",
    #             "/media/eze/Data/data/organismos/GCF_001624625.1/annotation/GCF_001624625.1_ASM162462v1_genomic.gb",
    #             tax=774)

    # mdb.delete_seq_collection("cruzi")
    # from_TriTrypDB(name="cruzi", tax=353153,
    #                fasta="/media/eze/Data/data/organismos/cruzi/TriTrypDB-34_TcruziCLBrenerEsmeraldo-like_Genome.fasta",
    #                gff="/media/eze/Data/data/organismos/cruzi/TriTrypDB-34_TcruziCLBrenerEsmeraldo-like.gff"
    #                )
    #_common_annotations("GCF_001624625.1","/tmp/GCF_001624625.1/")
    _common_annotations("cruzi","/tmp/cruzi/")
