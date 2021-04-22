"""

"""
import datetime
import logging as logging
import multiprocessing
import os
import pickle
import re
import shutil
import subprocess
from collections import defaultdict

import Bio.SearchIO as bpsio
import Bio.SeqIO as bpio
from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature, CompoundLocation
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from SNDG import init_log, mkdir, execute
from SNDG.BioMongo.Model.Feature import Feature, Location
from SNDG.BioMongo.Model.Protein import Protein
from SNDG.BioMongo.Model.SeqCollection import SeqCollection, AnnotationPipelineResult
from SNDG.BioMongo.Model.Sequence import BioProperty
from SNDG.BioMongo.Process.BioDocFactory import BioDocFactory
from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
from SNDG.BioMongo.Process.PathwaysAnnotator import PathwaysAnnotator
from SNDG.BioMongo.Process.SearchLoader import SearchLoader
from SNDG.BioMongo.Process.SearchLoader import load_hmm, load_blast_pdb

has_tax = False
try :
    from SNDG.BioMongo.Process.Taxon import Tax
    from SNDG.Sequence.ProteinAnnotator import Mapping
    has_tax = True
except:
    pass

from SNDG.Sequence import read_blast_table
from SNDG.Sequence import smart_parse as sp
from SNDG.Sequence.Hmmer import Hmmer

from SNDG.Sequence.so import SO_TERMS
from SNDG.WebServices.NCBI import NCBI
from SNDG.WebServices.Uniprot import Uniprot
from SNDG.Network.KEGG import Kegg

from SNDG.BioMongo.Process.Index import build_statistics, index_seq_collection

_log = logging.getLogger("Importer")


def _protein_iter(contigIterator, accept_feature=lambda f: ((f.type == "CDS)" and ("translation" in f.qualifiers))),
                  extract_annotation_feature=lambda f: f.type == "CDS",
                  extract_sequence=lambda c, f: f.qualifiers["translation"][
                      0] if "translation" in f.qualifiers else f.extract(c).seq.translate()):
    for c in contigIterator:
        for feature in c.features:
            if accept_feature(feature):
                f = extract_annotation_feature(feature)

                seq = extract_sequence(c, f)
                if "locus_tag" in f.qualifiers:
                    record = (SeqRecord(id=f.qualifiers["locus_tag"][0], description="",
                                        annotations=f.qualifiers, seq=seq), f)
                    if hasattr(feature, "sub_features"):
                        for sub_feature in feature.sub_features:
                            if sub_feature.type in ["transmembrane", "lipoprotein_signal_peptide", "signal_peptide"]:
                                start = (sub_feature.location.start - feature.location.start) / 3
                                stop = (sub_feature.location.end - feature.location.start) / 3
                                if feature.location.strand == -1:
                                    fstop = len(seq) - start
                                    fstart = len(seq) - stop
                                    sub_feature.location = FeatureLocation(fstart, fstop)
                                else:
                                    sub_feature.location = FeatureLocation(start, stop)
                                record[0].features.append(sub_feature)

                    yield record
                else:
                    _log.warning(str(f) + " has no locus_tag")


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


def from_ref_seq(name, ann_path, seqs=None, tax=None, tmp_dir=None,
                 extract_annotation_feature=lambda feature: feature.sub_features[
                     0] if feature.type == "gene" and hasattr(feature, "sub_features") and len(
                     feature.sub_features) else feature,
                 accept_protein_feature=lambda f: ((f.type == "CDS") and ("translation" in f.qualifiers)),
                 extract_sequence=lambda c, f: f.qualifiers["translation"][
                     0] if "translation" in f.qualifiers else f.extract(c).seq.translate(),
                 cpus=1):
    if seqs:
        seqs = {r.id: r.seq for r in bpio.parse(seqs, "fasta")}

    iter_seqs = list(sp(ann_path, seqs=seqs) if seqs else sp(ann_path))
    for contig in iter_seqs:
        if has_tax:
            seqCol = BioDocFactory.create_genome(name, contig, tax, Tax)
        else:
            seqCol = BioDocFactory.create_genome(name, contig)
        seqCol.save()
        break
    if not tmp_dir:
        tmp_dir = "/tmp/" + name + "/"
    mkdir(tmp_dir)
    gene_ids = {}
    with tqdm(iter_seqs) as pbar:
        for contig in pbar:
            pbar.set_description(contig.id)
            if len(contig.seq) > 15000000:
                contig.seq = ""
            contigDoc, gene_ids2 = BioDocFactory.create_contig(
                contig, seqCol, type_map=
                {"rRNA": "rRNA", "ncRNA": "ncRNA", NCBI.f_mRNA: NCBI.f_mRNA, "gene": "gene",
                 NCBI.f_CDS: NCBI.f_CDS, "rRNA": "rRNA", "tRNA": "tRNA", "tmRNA": "tmRNA"},
                extract_annotation_feature=extract_annotation_feature,
            )
            gene_ids.update(gene_ids2)
            contigDoc.save()

    prots = []

    with tqdm(_protein_iter(iter_seqs, accept_feature=accept_protein_feature,
                            extract_annotation_feature=extract_annotation_feature,
                            extract_sequence=extract_sequence
                            )) as pbar:
        for (protein, cds_f) in pbar:
            if "locus_tag" in cds_f.qualifiers:
                protDoc = BioDocFactory.create_protein(protein, cds_f)
                if len(protDoc.seq) > 30000:
                    raise Exception("No existen proteinas tan largas...")
                if protDoc.seq.count("*") > 1:
                    print (f"{cds_f.qualifiers['locus_tag'][0]}: Too many stop codons!")
                    continue
                if protDoc.seq.count("+") > 1:
                    print (f"{cds_f.qualifiers['locus_tag'][0]}: + signs found...!")
                    continue
                protDoc.gene_id = gene_ids[cds_f.qualifiers["locus_tag"][0]]
                protDoc.organism = name
                protDoc.auth = str(BioMongoDB.demo_id)
                protDoc.seq_collection_id = seqCol
                for f in protein.features:
                    protDoc.features.append(Feature(identifier=f.qualifiers["Ontology_term"][0], type=f.type,
                                                    location=Location(start=int(f.location.start),
                                                                      end=int(f.location.end))))

                prots.append(protDoc)
                if pbar.n and ((pbar.n % 1000) == 0):
                    Protein.objects.insert(prots)
                    prots = []
    if prots:
        Protein.objects.insert(prots)

    # _common_annotations(name, tmp_dir, cpu=cpus)
    return seqCol


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
            protDoc.auth = str(BioMongoDB.demo_id)
            protDoc.seq_collection_id = seqCol
            prots.append(protDoc)
            if pbar.n and ((pbar.n % 1000) == 0):
                Protein.objects.insert(prots)
                prots = []
    if prots:
        Protein.objects.insert(prots)

    _common_annotations(name, tmp_dir)


def create_proteome(tmp_dir, collection_name):
    protein_fasta = tmp_dir + "/proteins.fasta"
    if not os.path.exists(protein_fasta) or (not os.path.getsize(protein_fasta)):
        with open(protein_fasta, "w") as h:
            for p in Protein.objects(organism=collection_name).no_cache():
                bpio.write(SeqRecord(id=p.gene[0], description="", seq=Seq(p.seq)), h, "fasta")
    return protein_fasta


def _common_annotations_cmd(tmp_dir, protein_fasta, cpu=1, process_hmm=True, process_pdb=True):
    blast_result = None
    hmm_result = None
    if process_pdb:
        blast_result = tmp_dir + "/pdb_blast.xml"
        pdbs_path = "/data/databases/pdb/pdb_seqres.txt"
        if not os.path.exists(blast_result) or os.path.getsize(blast_result) < 10:
            cmd = f"blastp -qcov_hsp_perc 80 -max_hsps 1 -evalue 1e-5 -query {protein_fasta} -db {pdbs_path} -num_threads {cpu} -outfmt 5 -out {blast_result}"
            subprocess.call(cmd , shell=True)

    if process_hmm:
        hmm_result = tmp_dir + "/domains.hmm"
        params = {"--acc": None, "--cut_tc": None, "--notextw": None, "--cpu": str(cpu)}
        Hmmer(protein_fasta, output_file=hmm_result, params=params).query()

    return {"blast_pdb": blast_result, "hmm_result": hmm_result}


def common_annotations(collection_name, tmp_dir, cpu=1, remove_tmp=False):
    process_pdb = not Protein.objects(
        __raw__={"organism": collection_name, "features.type": SO_TERMS["polypeptide_structural_motif"]}).count()
    process_hmm = not (Protein.objects(__raw__={
        "organism": collection_name, "features.type": SO_TERMS["polypeptide_domain"]}).count())

    _common_annotations(collection_name, tmp_dir, cpu, remove_tmp, process_pdb, process_hmm)


def _common_annotations(collection_name, tmp_dir, cpu=1, remove_tmp=False, process_pdb=True, process_hmm=True):
    protein_fasta = create_proteome(tmp_dir, collection_name)

    results = _common_annotations_cmd(tmp_dir, protein_fasta, cpu, process_hmm, process_pdb)

    if process_pdb:
        blast_result = results["blast_pdb"]
        load_blast_pdb(collection_name, blast_result)
        if remove_tmp:
            if os.path.exists(blast_result):
                os.remove(blast_result)

    if process_hmm:
        hmm_result = results["hmm_result"]
        load_hmm(collection_name, hmm_result)
        if remove_tmp:
            if os.path.exists(hmm_result):
                os.remove(hmm_result)


def update_proteins(annotation_dir, proteome, seq_col_name, tax_id,
                    identity=0.9, cpus=multiprocessing.cpu_count(), db_init=None):

    # if db_init:
    #     from SNDG.Sequence.ProteinAnnotator import PABase
    #     PABase.sqldb.initialize(db_init)
    # mkdir(annotation_dir)
    # out = annotation_dir + "/species_blast.tbl"
    #
    # tax = Tax.select().where(Tax.ncbi_taxon_id == tax_id).get()
    # species_tax = None
    # for tax in Tax.parents(tax):
    #     if tax.node_rank == "genus":
    #         species_tax = tax
    #         break
    # tax_data = "/data/xomeq/tax/"
    # species_fasta = tax_data + str(int(species_tax.ncbi_taxon_id)) + ".fasta"

    if not os.path.exists(out):

        if not os.path.exists(species_fasta):
            Uniprot.download_proteome_from_tax(str(species_tax.ncbi_taxon_id), tax_data)

        cmd = "blastp -query %s  -db %s -evalue 0.00001 -outfmt 6  -max_hsps 1 -qcov_hsp_perc 0.9 -num_threads %i -out %s"
        execute(cmd % (
            proteome, species_fasta, cpus, out
        ))
    species_desc = {x.id.split("|")[1]: " ".join(x.description.split()[1:]) for x in bpio.parse(species_fasta, "fasta")}

    total = Protein.objects(organism=seq_col_name).count()
    with tqdm(bpsio.parse(out, "blast-tab"), total=total) as pbar:
        for query in pbar:
            pbar.set_description(query.id)
            if query[0][0].ident_pct > identity:

                unip = query[0].id.split("|")[1] if "|" in query[0].id else query[0].id
                dbxrefs = [x.db + "||" + x.value for x in Mapping.select().where(Mapping.uniprot == unip)]
                p = Protein.objects(gene=query.id, organism=seq_col_name).no_cache().get()

                if not p.description and unip in species_desc:
                    p.description = species_desc[unip].split("OS=")[0] + " | homology with: " + unip
                    p.save()

                if dbxrefs:
                    p = SearchLoader.update_protein_with_dbxref(
                        query.id, dbxrefs, seq_col_name)
                    p.save()


def load_pathways(genome_name, sbml_path, db, pathways_dir, prefered_biocyc=None,
                  gregexp="\(([\-\w\.]+)\)", gene_map=None, filter_file="allfilters_con_c.dat", sbml_desc=""):
    '''
    * gene_map: pikle file dict with sbml_gene_id as key and proteome_gene_id as value
    gregexp:  default = "\(([\-\w\.]+)\)"
    '''
    assert os.path.exists(sbml_path), sbml_path
    pathways_dir = pathways_dir + "/" if pathways_dir[-1] != "/" else pathways_dir
    if gene_map:
        assert os.path.exists(gene_map), "%s does not exists" % gene_map
        gene_map = pickle.load(open(gene_map))

    _log.info(db.proteins.update({"organism": genome_name, "properties._type": "pathways"},
                                 {"$pull": {"properties": {"_type": "pathways"}}}, multi=True))

    collection = SeqCollection.objects(name=genome_name).get()

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


def import_prop_blast(db, genome_name, offtarget_name, blast_output,
                      blast_output_type="table", description=None,
                      value_fn=lambda x: x, default_value=0, no_hit_value=0,
                      choices=[], type="number", defaultOperation=">"):
    genome = db.sequence_collection.find_one({"name": genome_name})
    drug_prop = [x for x in genome["druggabilityParams"] if x["name"] == offtarget_name]

    if drug_prop:
        drug_prop = drug_prop[0]
    else:
        if not description:
            description = offtarget_name
        drug_prop = {
            "target": "protein",
            "defaultGroupOperation": "max",
            "defaultValue": default_value,
            "name": offtarget_name,
            "defaultOperation": defaultOperation,
            "_cls": "SeqColDruggabilityParam",
            "uploader": "demo",
            "_class": "ar.com.bia.entity.SeqCollectionDoc",
            "type": type,
            "options": choices,
            "description": description
        }
        genome["druggabilityParams"].append(drug_prop)
        db.sequence_collection.save(genome)
    if blast_output_type == "xml":
        for query in bpsio.parse(blast_output, "blast-xml"):
            value = no_hit_value
            for hit in query:
                hsp = hit[0]
                value = value_fn(hsp)
            db.proteins.update({"organism": genome_name, "gene": query.id},
                               {"$set": {"search." + offtarget_name: value}})

    else:
        db.proteins.update({"organism": genome_name},
                           {"$set": {"search." + offtarget_name: no_hit_value}}, multi=True)
        for _, r in read_blast_table(blast_output).iterrows():
            db.proteins.update({"organism": genome_name, "gene": r["query"]},
                               {"$set": {"search." + offtarget_name: value_fn(r)}})




def import_kegg_annotation(db, genome_name, kegg_annotation):
    current_pathways = {}
    pw_kos = defaultdict(lambda: [])
    for ko_code, genes in tqdm(kegg_annotation.ko_gene.items()):
        if "ko:" + ko_code not in kegg_annotation.ko_dict:

            continue
        ko = kegg_annotation.ko_dict["ko:" + ko_code]
        if not ko:
            continue
        onts = ["kegg:" + ko_code]
        if ko["ecs"]:
            onts = onts + ko["ecs"]
        if ko["pathways"]:
            onts = onts + ko["pathways"]
            for x in ko["pathways"]:
                current_pathways[x] = 1
                pw_kos[x].append(ko_code)

        update = {
            "$addToSet": {"gene": {"$each": ko["genes"]},
                          "ontologies": {"$each": onts}}
        }
        # "$set": {"description": ko["desc"]} --> no  es precisa la desc, habria que ponerla solo en el caso que sea unknown function
        for g in genes:
            db.proteins.update({"organism": genome_name, "gene": g}, update, multi=True)

    keggs = []
    for pw in kegg_annotation.pw_dict.values():
        if pw["name"] in current_pathways:
            kos = list(set(pw_kos[pw["name"]]))
            kegg = {"term": pw["name"], "name": pw["title"], "properties": {"kos": kos}}
            if "reactions" in pw:
                kegg["count"] = len(pw["reactions"])
            keggs.append(kegg)

    db.sequence_collection.update({"name": genome_name}, {"$set": {"kegg": keggs}})


def fix_annotator_files(gff, fna, tag_replacement, tag="BIA_"):
    if not os.path.exists(fna + ".bk"):
        shutil.copy(fna, fna + ".bk")

    with open(fna, "w") as h:
        for x in bpio.parse(fna + ".bk", "fasta"):
            x.id = x.id.replace(tag, tag_replacement)
            x.name = x.id.replace(tag, tag_replacement)
            x.description = ""
            bpio.write(x, h, "fasta")

    contigs = []

    if not os.path.exists(gff + ".bk"):
        shutil.copy(gff, gff + ".bk")

    with open(gff + ".bk") as r:
        with open(gff, "w") as w:
            data = r.read().split("##FASTA")
            for i, gffstrs in enumerate(data[:-1] if len(data) > 1 else data):
                if i == 0:
                    w.write(gffstrs)
                else:
                    w.write(gffstrs.split("##feature-ontology so.obo")[1])
    data = list(GFF.parse(gff))
    with open(gff, "w") as h:
        for x in data:
            if not len([y for y in x.features if y.type == "contig"]):
                continue
            x.id = x.id.replace(tag, tag_replacement)
            x.description = ""
            for f in x.features:
                if "locus_tag" in f.qualifiers:
                    f.qualifiers["locus_tag"] = [x.id + "_" + f.qualifiers["locus_tag"][0].replace(tag, "")]
                for sf in f.sub_features:
                    if "locus_tag" in sf.qualifiers:
                        sf.qualifiers["locus_tag"] = [x.id + "_" + sf.qualifiers["locus_tag"][0].replace(tag, "")]
            contigs.append(x)
        GFF.write(contigs, h, False)


if __name__ == '__main__':
    init_log()

    logging.getLogger("peewee").setLevel(logging.WARN)
    from peewee import MySQLDatabase
    from SNDG.BioMongo.Process.Taxon import tax_db

    tax_db.initialize(MySQLDatabase('bioseqdb', user='root', passwd="mito"))
    mdb = BioMongoDB("saureus", port=27017)


    # mdb.delete_seq_collection("ILEX_PARA2")

    def extract_annotation_feature(feature):
        mrnas = [f for f in feature.sub_features if f.type == "mRNA"]
        return mrnas[
            0] if feature.type == "gene" and len(mrnas) else feature


    def accept_protein_feature(feature):
        return feature.type == "gene" and feature.sub_features and feature.sub_features[0].type == "mRNA"


    # prot_dict = bpio.to_dict(bpio.parse("/data/organismos/ILEX_PARA/contigs/ncbi_IP4.faa","fasta"))
    def extract_sequence(c, f):
        return prot_dict[f.id].seq
        # def getphase(pfeature):
        #     phase = 0
        #     if "phase" in pfeature.qualifiers:
        #         phase = int(pfeature.qualifiers["phase"][0])
        #     return phase
        # locs = sorted([(x.location,getphase(x))  for x in f.sub_features if x.type == "CDS"],key=lambda x:x[0].start)
        #
        # if len(locs) > 1:
        #     if f.strand == -1:
        #         locs =  list(reversed(locs))
        #     if locs[0][1] != 0:
        #
        #         locs[0][0] = FeatureLocation(start= locs[0][0].start + locs[0][1] if f.strand == 1 else  locs[0][0].start
        #                                      ,end= locs[0][0].end - locs[0][1] if f.strand == -1 else  locs[0][0].end
        #                                      ,strand=locs[0].strand)
        #     l = CompoundLocation([x[0] for x in locs])
        # else:
        #     l = locs[0][0]
        # cds = SeqFeature(strand=f.strand, location=l)
        # seq = cds.extract(c).seq.translate()
        # if len(seq) == 0:
        #     log.warning("empty sequence!")
        #
        # if seq[-1] == "*":
        #     seq  = seq[:-1]
        # if "*" in str(seq):
        #     _log.warning("stop codons in the middle of the sequence")
        # return seq


    # from_ref_seq("ILEX_PARA2", "/data/organismos/ILEX_PARA/contigs/ncbi_IP4.gff3",
    #              extract_annotation_feature=extract_annotation_feature,
    #              accept_protein_feature=accept_protein_feature,
    #              extract_sequence=extract_sequence,
    #              seqs="/data/organismos/ILEX_PARA/contigs/ncbi_IP4.fna", tax=185542,
    #              tmp_dir="/data/organismos/ILEX_PARA/tmp", cpus=4)
    # _common_annotations("ILEX_PARA2", "/data/organismos/ILEX_PARA2/", cpu=1)
    # from SNDG.BioMongo.Process.Index import index_seq_collection,build_statistics
    # index_seq_collection(mdb.db, "ILEX_PARA2", keywords=True, pathways=False, structure=True)
    # build_statistics(mdb.db,"ILEX_PARA2")

    # mysqldb = ProteinAnnotator.connect_to_db(database="unipmap", user="root", password="mito")

    seq_col_name = "ILEX_PARA2"
    # kegg_annotation = Kegg()
    # # kegg_annotation.update_files()
    # kegg_annotation.init()
    # ilex_data = "/data/organismos/ILEX_PARA2/query.ko"
    # kegg_annotation.read_annotation(ilex_data)
    # import_kegg_annotation(mdb.db,seq_col_name, kegg_annotation)

    # from SNDG.BioMongo.Process.PathwaysAnnotator import PathwaysAnnotator
    # import pymongo
    # pa = PathwaysAnnotator(mdb.db, "ILEX_PARA_TRANSCRIPT", "/data/organismos/ILEX_PARA_TRANSCRIPT/sbml/")
    # pa.sbml("/data/organismos/ILEX_PARA_TRANSCRIPT/sbml/ilex_para.sbml")
    # pa.annotate()

    filter_tax = {2,
                  22,
                  29,
                  81,
                  192,
                  193,
                  194,
                  286,
                  356,
                  396,
                  429,
                  434,
                  482,
                  506,
                  543,
                  629,
                  641,
                  666,
                  724,
                  810,
                  816,
                  817,
                  838,
                  914,
                  919,
                  976,
                  1046,
                  1090,
                  1117,
                  1118,
                  1224,
                  1227,
                  1236,
                  1239,
                  1279,
                  1297,
                  1760,
                  1762,
                  1763,
                  1883,
                  2037,
                  2062,
                  2063,
                  2093,
                  28211,
                  28216,
                  28221,
                  32011,
                  32066,
                  33877,
                  33882,
                  35798,
                  39773,
                  40117,
                  40222,
                  44249,
                  55158,
                  67814,
                  68297,
                  68336,
                  72276,
                  74152,
                  82115,
                  85006,
                  85007,
                  85008,
                  85011,
                  85023,
                  85025,
                  91061,
                  91347,
                  93681,
                  135613,
                  135617,
                  135623,
                  142182,
                  171552,
                  186801,
                  186826,
                  188708,
                  191411,
                  191412,
                  200783,
                  200795,
                  200918,
                  200930,
                  200938,
                  200940,
                  201174,
                  203682,
                  203691,
                  204441,
                  204455,
                  204457,
                  267890,
                  508458,
                  544448,
                  578822,
                  1293497,
                  2157}

    # from SNDG.BioMongo.Process.BioCyc2Mongo import BioCyc
    # biocyc = BioCyc("saureus")
    # biocyc.complete_pathways("ILEX_PARA_TRANSCRIPT", "/data/databases/biocyc/metacyc/pathways.dat",
    #                           "/data/databases/biocyc/metacyc/reactions.dat", filter_tax)

    # index_seq_collection(mdb.db, "ILEX_PARA_TRANSCRIPT",ec=False, go=True, keywords=True, organism_idx=True, pathways=True,
    #                      structure=True)
    # build_statistics(mdb.db,"ILEX_PARA_TRANSCRIPT")

    # base = "/data/projects/ecoli_mex/"

    # for x in tqdm(["AIEC_C4435", "AIEC_C7225", "AIEC_C7230", "AIEC_C7223",
    #                "AIEC_C7226"]):  # ,"160A" "188B","1CT136A", "22A", "86A"
    #     mdb.protein_fasta("/home/eze/Downloads/" + x + ".fasta", "Eco" + x)
    #     kegg_annotation = Kegg()
    #     kegg_annotation.init()
    #     ilex_data = "/data/organismos/Eco" + x + "/annotation/query.ko"
    #     kegg_annotation.read_annotation(ilex_data)
    #     import_kegg_annotation(mdb.db, "Eco" + x, kegg_annotation)

    # mdb.delete_seq_collection("Eco" + x)
    # gff = glob.glob(base + x + "/ncbi_*.gff3")[0]
    # fna = glob.glob(base + x + "/ncbi_*.fna")[0]
    #
    # fix_annotator_files(gff,fna,"Eco" + x)
    #
    #
    # def extract_annotation_feature(feature):
    #     mrnas = [f for f in feature.sub_features if f.type == "mRNA"]
    #     return mrnas[
    #         0] if feature.type == "gene" and len(mrnas) else feature
    #
    #
    # seqCol = from_ref_seq("Eco" + x, base + x + ".gbk", tax=562, tmp_dir="/tmp/Eco" + x + "", cpus=3)
    #
    # seqCol.description = "Escherichia coli " + x + ""
    # seqCol.organism = "Escherichia coli " + x + ""
    # seqCol.auth = "40"
    # seqCol.save()
    # mdb.protein_fasta("/tmp/Eco" + x + "/genome.fasta", "Eco" + x)
    # update_proteins("/tmp/Eco" + x + "/", "/tmp/Eco" + x + "/genome.fasta", "Eco" + x, 562, db_init=mysqldb)
    # from SNDG.WebServices.Offtargeting import Offtargeting
    #
    # Offtargeting.offtargets("/tmp/Eco" + x + "/genome.fasta",
    #                         "/data/organismos/Eco" + x + "/annotation/offtarget/"
    #                         )
    #
    # import_prop_blast(mdb.db, "Eco" + x, "hit_in_deg",
    #                   "/data/organismos/Eco" + x + "/annotation/offtarget/degaa-p.tbl",
    #                   "table", "Hit in DEG database",
    #                   value_fn=lambda x: x.identity > 70,
    #                   default_value=True,
    #                   no_hit_value=False, choices=["True", "False"], type="value", defaultOperation="equal")
    #
    # import_prop_blast(mdb.db, "Eco" + x, "human_offtarget",
    #                   "/data/organismos/Eco" + x + "/annotation/offtarget/gencode.tbl",
    #                   "table", "Human offtarget score (1 - best hit identity)",
    #                   value_fn=lambda x: 1 - (x.identity * 1.0 / 100),
    #                   default_value=0.4,
    #                   no_hit_value=1)
    # import_prop_blast(mdb.db, "Eco" + x, "gut_microbiota_offtarget",
    #                   "/data/organismos/Eco" + x + "/annotation/offtarget/gut_microbiota.tbl",
    #                   "table", "Gut microbiota offtarget score (1 - best hit identity)",
    #                   value_fn=lambda x: 1 - (x.identity * 1.0 / 100),
    #                   default_value=0.4,
    #                   no_hit_value=1)

    # index_seq_collection(mdb.db,"Eco" + x,pathways=False,go=True,keywords=True,ec=True,organism_idx=True,structure=True)
    # build_statistics(mdb.db,"Eco" + x)

    # import logging
    #

    # from SNDG.BioMongo.Process.Index import index_seq_collection
    #

    #
    # tofix = [u'19', u'23', u'36', u'43', u'54', u'64', u'GCF_000508085.1', u'GCF_000373365.1', u'GCF_000966285.1',
    #          u'GCF_000805695.1', u'GCF_001580035.1', u'GCF_000333795.1', u'GCF_000788295.1', u'GCF_000510935.1',
    #          u'GCF_000769675.1', u'GCF_000188155.2', u'GCF_000213355.1', u'GCF_000179595.2', u'GCF_000213335.1',
    #          u'GCF_000213395.1', u'GCF_002013775.1', u'GCF_002013745.1', u'GCF_002013685.1', u'GCF_002013645.1',
    #          u'GCF_002013545.1']
    # for seq_col_name in tqdm(tofix):
    #     tid = int(mdb.db.sequence_collection.find_one({"name": seq_col_name})["tax"]["tid"])
    #     tmp_dir = "/data/organismos/" + seq_col_name + "/annotation/"
    #     proteome_dir = "/data/organismos/" + seq_col_name + "/contigs/"
    #     mkdir(tmp_dir)
    #     mkdir(proteome_dir)
    #     protein_fasta = create_proteome(proteome_dir, seq_col_name)
    #     update_proteins(tmp_dir, protein_fasta, seq_col_name, tid)
    #
