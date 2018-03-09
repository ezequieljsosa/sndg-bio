"""

"""
import logging as logging
import os
import pickle
import re
import subprocess
from collections import defaultdict

import Bio.SeqIO as bpio
import datetime
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from SNDG import init_log
from SNDG import mkdir
from SNDG.BioMongo.Model.Protein import Protein
from SNDG.BioMongo.Model.SeqCollection import SeqCollection, AnnotationPipelineResult
from SNDG.BioMongo.Model.Sequence import BioProperty
from SNDG.BioMongo.Process.BioDocFactory import BioDocFactory
from SNDG.BioMongo.Process.PathwaysAnnotator import PathwaysAnnotator
from SNDG.BioMongo.Process.SearchLoader import load_hmm, load_blast_pdb
from SNDG.BioMongo.Process.Taxon import Tax
from SNDG.Sequence import smart_parse as sp
from SNDG.Sequence.Hmmer import Hmmer
from SNDG.Sequence.so import SO_TERMS
from SNDG.WebServices.NCBI import NCBI

_log = logging.getLogger("Importer")


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
                        seq = cds.extract(c).seq.translate()
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

    _common_annotations(name, tmp_dir)


def from_TriTrypDB(name, gff, fasta, tax, tmp_dir=None):
    genome = {x.id: x for x in sp(fasta)}
    from BCBio import GFF
    import re
    annotation = list(GFF.parse(gff, base_dict=genome))
    contig = annotation[0]

    seqCol = BioDocFactory.create_genome(name, contig, tax, Tax)
    seqCol.save()

    if not tmp_dir:
        tmp_dir = "/tmp/" + name + "/"
    mkdir(tmp_dir)
    gene_ids = {}
    with tqdm(annotation) as pbar:
        for contig in pbar:
            pbar.set_description(contig.id)
            if len(contig.seq) > 15000000:
                contig.seq = ""
            contigDoc, gene_ids2 = BioDocFactory.create_contig(contig, seqCol, type_map=
            {"rRNA": "rRNA", "ncRNA": "ncRNA", NCBI.f_mRNA: "gene", "exon": "exon", "gene": "gene",
             NCBI.f_CDS: NCBI.f_CDS, "rRNA": "rRNA", "tRNA": "tRNA", "tmRNA": "tmRNA",
             "snoRNA": "snoRNA",
             "three_prime_UTR": "three_prime_UTR", "five_prime_UTR": "five_prime_UTR"
             })
            gene_ids.update(gene_ids2)
            contigDoc.save()
    prots = []
    with tqdm(tritryp_protein_iter(annotation)) as pbar:
        for (protein, cds_f) in pbar:

            protDoc = Protein(seq=str(protein.seq), name=protein.id)

            if "description" in cds_f.qualifiers:
                protein_description = cds_f.qualifiers['description'][0]
            elif "Note" in cds_f.qualifiers:
                protein_description = cds_f.qualifiers['Note'][0]
            elif "product" in cds_f.qualifiers:
                protein_description = cds_f.qualifiers['product'][0]
            else:
                protein_description = ""

            protDoc.description = protein_description

            gos = []
            if "Ontology_term" in cds_f.qualifiers:
                gos = [x.lower() for x in cds_f.qualifiers["Ontology_term"] if
                       "GO:" in x and (x not in ["GO:0008150", "GO:0003674", "GO:0005575"])]

            note = cds_f.qualifiers["Note"][0].split(" ")[0] if "Note" in cds_f.qualifiers else ""
            ecs = ["ec:" + note] if re.match('^[0-9]+\.[0-9\-]+\.[0-9\-]+\.[0-9\-]$', note) else []
            ontologies = list(set(ecs + gos))

            protDoc.gene = [protein.id]
            protDoc.ontologies = ontologies
            protDoc.alias = [protein.id]

            if len(protDoc.seq) > 30000:
                raise Exception("No existen proteinas tan largas...")
            protDoc.gene_id = gene_ids[protein.id]
            protDoc.organism = name
            protDoc.auth = BioMongoDB.demo_id
            protDoc.seq_collection_id = seqCol
            prots.append(protDoc)
            if pbar.n and ((pbar.n % 1000) == 0):
                Protein.objects.insert(prots)
                prots = []
    if prots:
        Protein.objects.insert(prots)

    _common_annotations(name, tmp_dir)


def _common_annotations(name, tmp_dir):
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
        # os.remove(blast_result)

    if not Protein.objects(__raw__={"organism": name, "features.type": SO_TERMS["polypeptide_domain"]}).count():
        hmm_result = tmp_dir + "/domains.hmm"
        Hmmer(protein_fasta, output_file=hmm_result).query()
        load_hmm(name, hmm_result)
        os.remove(hmm_result)


def load_pathways( genome_name, sbml_path,db,pathways_dir, prefered_biocyc=None,
                  gregexp="\(([\-\w\.]+)\)", gene_map=None, filter_file="allfilters_con_c.dat", sbml_desc=""):
    '''
    * gene_map: pikle file dict with sbml_gene_id as key and proteome_gene_id as value
    gregexp:  default = "\(([\-\w\.]+)\)"
    '''
    assert os.path.exists(sbml_path), sbml_path
    if gene_map:
        assert os.path.exists(gene_map), "%s does not exists" % gene_map
        gene_map = pickle.load(open(gene_map))

    _log.info(db.proteins.update({"organism": genome_name, "properties._type": "pathways"},
                                      {"$pull": {"properties": {"_type": "pathways"}}}, multi=True))

    collection = SeqCollection.objects(name=genome_name).get()

    _log.info(db.proteins.update({"organism": genome_name, "properties._type": "pathways"},
                                      {"$pull": {"properties": {"_type": "pathways"}}}, multi=True))
    _log.info(
        db.proteins.update({"organism": genome_name, "reactions": {"$exists": True}}, {"$set": {"reactions": []}},
                                multi=True))

    ps = PathwaysAnnotator(db, genome_name, pathways_dir)
    ps.sbml(sbml_path)
    ps.species_filter(filter_file)
    gene_name_regexp = re.compile(gregexp)
    if gene_map:
        def map_gene(notes):
            genes = gene_name_regexp.findall(notes)
            return [gene_map[x] if x in gene_map else x for x in genes]

        ps.extract_genes_from_notes(map_gene)
    else:
        ps.extract_genes_from_notes(lambda notes: gene_name_regexp.findall(notes))
    if prefered_biocyc:
        ps.prefered_biocyc = prefered_biocyc
    ps.annotate()

    result = AnnotationPipelineResult(name="pathways", version=1, timestamp=datetime.datetime.now)
    result.inputs = {
        "prefered_biocyc": prefered_biocyc,
        "gregexp": gregexp,
        "gene_map": gene_map,
        "sbml": gene_map,
        "sbml_desc": sbml_desc
    }

    chokepoints = db.proteins.count({"organism": genome_name, "properties.property": "chokepoint"})

    result.results = {
        "pathways": len(ps.sbmlprocessor.pathways),
        "annotated_proteins": db.proteins.count({"organism": genome_name, "reactions.0": {"$exists": True}}),
        "chokepoints": chokepoints,
        "unmapped_genes": ps.unmapped_genes
    }

    _log.info("PW annotation finished")

    collection.pathways = ps.pathways
    collection.pipelines.append(result)

    collection.save()


def correct_chokes(self, name):
    metabolites_in = defaultdict(list)
    metabolites_out = defaultdict(list)
    for p in self.db.proteins.find({"organism": name, "reactions.0": {"$exists": True}}):
        for r in p["reactions"]:
            for m in r["products"]:
                metabolites_out[m["name"]].append(r["name"])
            for m in r["substrates"]:
                metabolites_in[m["name"]].append(r["name"])

    for m, r in metabolites_in.items():
        if (len(set(r)) > 1):
            # or (self.db.proteins.count({"organism":name,"reactions.name": r[0]}) > 1)):
            del metabolites_in[m]
    for m, r in metabolites_out.items():
        if (len(set(r)) > 1):
            # or (self.db.proteins.count({"organism":name,"reactions.name": r[0]}) > 1)):
            del metabolites_out[m]

    choke_reactions_in = []
    for rs in metabolites_in.values():
        choke_reactions_in += rs
    choke_reactions_out = []
    for rs in metabolites_out.values():
        choke_reactions_out += rs

    reaction_metabolites = defaultdict(lambda: [])
    for m, rs in metabolites_in.items():
        for r in rs:
            reaction_metabolites[r].append(m)
    for m, rs in metabolites_out.items():
        for r in rs:
            reaction_metabolites[r].append(m)

    for p in Protein.objects(organism=name, reactions__0__exists=True).no_cache().timeout(False):
        cout = bool([r.name for r in p.reactions if r.name in choke_reactions_out])
        cin = bool([r.name for r in p.reactions if r.name in choke_reactions_in])
        p.search.chokepoint = cout | cin
        if p.search.chokepoint:
            p.search.chokepoint_type = "double" if (cout & cin) else ("production" if cout else "consuming")
            prop = [x for x in p.properties if x.property == "chokepoint"]
            if prop:
                prop = prop[0]
                prop.metabolites = []
                for x in p.reactions:
                    if x.name in reaction_metabolites:
                        prop.metabolites += reaction_metabolites[x.name]
            else:
                prop = BioProperty(_type="pathways", property="chokepoint", metabolites=[],
                                   type=p.search.chokepoint_type)
                for x in p.reactions:
                    if x.name in reaction_metabolites:
                        prop.metabolites += reaction_metabolites[x.name]
                p.properties.append(prop)
        else:
            del p.search.chokepoint_type
            p.properties = [x for x in p.properties if x.property != "chokepoint"]
        p.save()


if __name__ == '__main__':
    init_log()
    from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
    import logging

    logging.getLogger("peewee").setLevel(logging.WARN)
    from peewee import MySQLDatabase
    from SNDG.BioMongo.Process.Taxon import tax_db

    tax_db.initialize(MySQLDatabase('bioseqdb', user='root', passwd="mito"))
    mdb = BioMongoDB("tdr")
    mdb.delete_seq_collection("GCF_001624625.1")
    from_ref_seq("GCF_001624625.1",
                 "/data/organismos/GCF_001624625.1/annotation/GCF_001624625.1_ASM162462v1_genomic.gb",
                 tax=774)



    mdb.delete_seq_collection("cruzi")
    from_TriTrypDB(name="cruzi", tax=353153,
                   fasta="/data/organismos/cruzi/TriTrypDB-34_TcruziCLBrenerEsmeraldo-like_Genome.fasta",
                   gff="/data/organismos/cruzi/TriTrypDB-34_TcruziCLBrenerEsmeraldo-like.gff"
                   )

    # _common_annotations("cruzi","/tmp/cruzi/")
    load_pathways("cruzi", "/data/organismos/cruzi/pathways/pathways-sm.sbml",
                  pymongo.MongoClient().tdr,"/data/organismos/cruzi/pathways/")
    load_pathways("GCF_001624625.1", "/data/organismos/GCF_001624625.1/pathways/pathways-sm.sbml",
              pymongo.MongoClient().tdr,"/data/organismos/GCF_001624625.1/pathways/")
    load_pathways("Pext14-3B", "/data/organismos/Pext14-3B/pathways/pathways-sm.sbml",
                  pymongo.MongoClient().tdr,"/data/organismos/Pext14-3B/pathways/")

