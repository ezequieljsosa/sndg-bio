"""

"""
import logging as logging
import os
import multiprocessing
import pickle
import re
import subprocess
from collections import defaultdict

import Bio.SeqIO as bpio
import Bio.SearchIO as bpsio
import datetime
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from SNDG import init_log ,mkdir ,execute

from SNDG.BioMongo.Model.Protein import Protein
from SNDG.BioMongo.Model.SeqCollection import SeqCollection, AnnotationPipelineResult
from SNDG.BioMongo.Model.Sequence import BioProperty
from SNDG.BioMongo.Process.BioDocFactory import BioDocFactory
from SNDG.BioMongo.Process.PathwaysAnnotator import PathwaysAnnotator
from SNDG.BioMongo.Process.SearchLoader import load_hmm, load_blast_pdb
from SNDG.BioMongo.Process.Taxon import Tax
from SNDG.Sequence import smart_parse as sp
from SNDG.Sequence.Hmmer import Hmmer
from SNDG.Sequence.ProteinAnnotator import ProteinAnnotator, Mapping
from SNDG.Sequence.so import SO_TERMS
from SNDG.WebServices.NCBI import NCBI
from SNDG.WebServices.Uniprot import Uniprot
from SNDG.BioMongo.Process.SearchLoader import SearchLoader


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

def create_proteome(tmp_dir,collection_name):
    protein_fasta = tmp_dir + "/proteins.fasta"
    if not os.path.exists(protein_fasta) or (not os.path.getsize(protein_fasta)):
        with open(protein_fasta, "w") as h:
            for p in Protein.objects(organism=collection_name).no_cache():
                bpio.write(SeqRecord(id=p.gene[0], description="", seq=Seq(p.seq)), h, "fasta")
    return protein_fasta


def _common_annotations_cmd(name, tmp_dir,protein_fasta, cpu=1,process_hmm=True,process_pdb=True):
    if process_pdb:
        blast_result = tmp_dir + "/pdb_blast.xml"
        pdbs_path = "/data/databases/pdb/processed/seqs_from_pdb.fasta"
        if not os.path.exists(blast_result):
            cmd = "blastp -qcov_hsp_perc 80 -max_hsps 1 -evalue 1e-5 -query %s -db %s -num_threads %i -outfmt 5 -out %s"
            subprocess.call(cmd % (protein_fasta, pdbs_path, cpu, blast_result), shell=True)

    if process_hmm:
        hmm_result = tmp_dir + "/domains.hmm"
        params = {"--acc": None, "--cut_tc": None, "--notextw": None, "--cpu": cpu}
        Hmmer(protein_fasta, output_file=hmm_result,params=params).query()


def _common_annotations(collection_name, tmp_dir, cpu=1, remove_tmp=False):
    process_pdb = not Protein.objects(
        __raw__={"organism": collection_name, "features.type": SO_TERMS["polypeptide_structural_motif"]}).count()
    process_hmm = not Protein.objects(__raw__={
        "organism": collection_name, "features.type": SO_TERMS["polypeptide_domain"]}).count()

    protein_fasta= create_proteome(tmp_dir,collection_name)

    _common_annotations_cmd(collection_name, tmp_dir,protein_fasta, cpu,process_hmm,process_pdb)


    if process_pdb:
        load_blast_pdb(collection_name, blast_result)
        if remove_tmp:
            if os.path.exists(blast_result):
                os.remove(blast_result)

    if process_hmm:
        load_hmm(collection_name, hmm_result)
        if remove_tmp:
            if os.path.exists(hmm_result):
                os.remove(hmm_result)



def update_proteins(annotation_dir, proteome,seq_col_name, tax_id, cpus=multiprocessing.cpu_count()):

    mkdir(annotation_dir)
    out = annotation_dir + "/species_blast.tbl"

    if not os.path.exists(out):
        tax = Tax.select().where(Tax.ncbi_taxon_id == tax_id).get()
        species_tax = None
        for tax in Tax.parents(tax):
            if tax.node_rank == "genus":
                species_tax = tax
                break

        tax_data = "/data/xomeq/tax/"
        species_fasta =tax_data  + str(int(species_tax.ncbi_taxon_id)) + ".fasta"
        if not os.path.exists(species_fasta):
            Uniprot.download_proteome_from_tax(str(species_tax.ncbi_taxon_id), tax_data)


        cmd = "blastp -query %s  -db %s -evalue 0.00001 -outfmt 6  -max_hsps 1 -qcov_hsp_perc 0.9 -num_threads %i -out %s"
        execute(cmd % (
            proteome, species_fasta, cpus, out
        ))

    ProteinAnnotator.connect_to_db(database="unipmap", user="root", password="mito")

    with tqdm(list(bpsio.parse(out, "blast-tab"))) as pbar:
        for query in pbar:
            pbar.set_description(query.id)
            try:
                if query[0][0].ident_pct > 0.9:
                    unip = query[0].id.split("|")[1] if "|" in query[0].id else query[0].id
                    dbxrefs = [x.db + "||" + x.value for x in Mapping.select().where(Mapping.uniprot == unip)]
                    p = SearchLoader.update_protein_with_dbxref(query.id, dbxrefs, seq_col_name)
                    p.save()
            except Exception as e:
                _log.warn(e)


def load_pathways(genome_name, sbml_path, db, pathways_dir, prefered_biocyc=None,
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
    import pymongo

    logging.getLogger("peewee").setLevel(logging.WARN)
    from peewee import MySQLDatabase
    from SNDG.BioMongo.Process.Taxon import tax_db
    from SNDG.BioMongo.Process.Index import index_seq_collection,build_statistics

    #tax_db.initialize(MySQLDatabase('bioseqdb', user='root', passwd="mito"))
    mdb = BioMongoDB("saureus")

    tofix = [u'19', u'23', u'36', u'43', u'54', u'64', u'GCF_000508085.1', u'GCF_000373365.1', u'GCF_000966285.1',
             u'GCF_000805695.1', u'GCF_001580035.1', u'GCF_000333795.1', u'GCF_000788295.1', u'GCF_000510935.1',
             u'GCF_000769675.1', u'GCF_000188155.2', u'GCF_000213355.1', u'GCF_000179595.2', u'GCF_000213335.1',
             u'GCF_000213395.1', u'GCF_002013775.1', u'GCF_002013745.1', u'GCF_002013685.1', u'GCF_002013645.1',
             u'GCF_002013545.1']
    for seq_col_name in tdqm(tofix):

        tid = int(mdb.db.sequence_collection.find_one({"name":seq_col_name})["tax"]["tid"])
        tmp_dir = "/data/organismos/" + seq_col_name + "/annotation/"
        proteome_dir = "/data/organismos/" + seq_col_name + "/contigs/"
        mkdir(tmp_dir)
        mkdir(proteome_dir)
        protein_fasta= create_proteome(proteome_dir,seq_col_name)
        update_proteins(tmp_dir, protein_fasta,seq_col_name, tid )

        index_seq_collection(mdb.db,seq_col_name,pathways=False,structure=False)
        build_statistics(mdb.db,seq_col_name)