"""

"""

import logging
import os
import re

from bson.objectid import ObjectId
from tqdm import tqdm

import Bio.SearchIO as bpsio
from SNDG.WebServices.Offtarget import Offtarget
from SNDG import mkdir
from SNDG.BioMongo.Model.Alignment import SimpleAlignment, AlnLine
from SNDG.BioMongo.Model.Feature import Feature, Location
from SNDG.BioMongo.Model.Protein import Protein
from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
from SNDG.BioMongo.Process.BioDocFactory import BioDocFactory
from SNDG.Sequence import identity, coverage, hit_coverage
from SNDG.Sequence.so import SO_TERMS

_log = logging.getLogger(__name__)


class SearchLoader():

    @staticmethod
    def offtarget(organism,offtarget_databases,offtarget_names,tmp_dir=None):
        if not tmp_dir:
            tmp_dir = "/data/organismos/" + organism + "/annotation/"
        mkdir(tmp_dir)
        proteins = tmp_dir + "proteins.fasta"
        if not os.path.exists(proteins):
            BioMongoDB.protein_fasta(proteins,organism)
        results = Offtarget.offtargets(proteins,tmp_dir,offtarget_databases)
        for i,name in enumerate(offtarget_names):
            load_blast_features( organism, results[i], name,
                                 min_identity=0.4, min_query_coverage=0.4, min_hit_coverage=0.4)




    @staticmethod
    def update_protein_with_dbxref(protein_gene, annotations,organism):
        p = Protein.objects(gene=protein_gene,organism=organism).no_cache().get()

        for ann in annotations:
            word = ann.split("||")[1]
            if ann.lower().startswith("ec") or ann.lower().startswith("go"):
                p.ontologies.append(word)
                p.keywords.append(word)
            for attr in """EcoGene
            Ensembl_PRO
            Ensembl_TRS
            WormBase_PRO
            WormBase_TRS
            UniGene
            GeneDB
            EuPathDB""".split():
                if ann.startswith(attr):
                    p.keywords.append(word)
                    p.alias.append(word)
            for g in ["Gene_Name","Gene_OrderedLocusName","Gene_ORFName"]:
                if ann.startswith(g):
                    p.keywords.append(word)
                    p.alias.append(word)
                    p.gene.append(word)

        #if "Uncharacterized" in p.description:
        #    p.description =
        p.ontologies = list( set([x.lower() for x in p.ontologies]))
        p.keywords = list(set([x.lower() for x in p.keywords]))
        p.alias = list(set(p.alias))
        p.gene = [protein_gene] +  [x for x in set(p.gene) if x != protein_gene]
        return p




def load_hmm(organism, hmm_file, transform_query_regexp=None, transform_hit_regexp=None):
    assert os.path.exists(hmm_file)
    for query in tqdm(bpsio.parse(hmm_file, 'hmmer3-text')):
        for hit in query:
            for hsp in hit:
                gene = query.id
                if transform_query_regexp:
                    gene = re.search(transform_query_regexp, query.id, re.IGNORECASE).group(1)

                hit_name = hit.id
                if transform_query_regexp:
                    hit_name = re.search(transform_hit_regexp, hit_name, re.IGNORECASE).group(1)

                proteins = Protein.objects(organism=organism,alias=gene).no_cache().timeout(False)
                for protein in proteins:
                    dn = [d for d in protein.domains() if
                          (d.identifier == hit_name) and (d.location.start == hsp.query_start) and (
                                  d.location.end == hsp.query_end)]
                    if dn:
                        protein.features.remove(dn[0])

                    hsp_feature = Feature(_id=ObjectId(), location=Location(start=hsp.query_start, end=hsp.query_end),
                                          aln=SimpleAlignment(
                                              evalue=hsp.evalue,
                                              aln_query=AlnLine(name=hsp.query_id, seq=str(hsp.aln[0].seq),
                                                                start=hsp.query_start, end=hsp.query_end),
                                              aln_hit=AlnLine(name=hsp.hit.id, seq=str(hsp.aln[1].seq),
                                                              start=hsp.hit_start, end=hsp.hit_end),
                                              aln_cd=hsp.aln_annotation["CS"] if "CS" in hsp.aln_annotation else "",
                                              aln_pp=hsp.aln_annotation["PP"] if "PP" in hsp.aln_annotation else ""
                                          ),
                                          identifier=hsp.hit.id, type=SO_TERMS["polypeptide_domain"])

                    protein.features.append(hsp_feature)
                    protein.save()


def load_blast_pdb(organism, blast_file):
    '''
    load_blast_features.py -db tdr -org TGONDII -ft SO:0001079 -i 0.95  -qc 0.8 -hc 0.8 -b ./blast_pdb.xml
    '''
    load_blast_features(organism, blast_file, SO_TERMS["polypeptide_structural_motif"], 0.95, 0.8, 0.8)


def load_pdb_domains(organism, blast_file, feature_type="SO:0001079",
                     min_identity=0.9, min_query_coverage=0.9, min_hit_coverage=0.9):
    queries = list(bpsio.parse(blast_file, 'blast-xml'))
    features_added = 0
    with tqdm(queries) as pbar:
        for query in queries:
            pbar.set_description(query.id)
            pfam, dnstart, dnend = query.id.split("_")[-3:]
            dnstart, dnend = int(dnstart), int(dnend)
            gene = "_".join(query.id.split("_")[:-3])

            proteins = Protein.objects(organism=organism, gene=gene).no_cache().timeout(False)
            change = False
            for protein in proteins:
                for hit in query:
                    hsp = hit[0]
                    dn = [x for x in protein.domains() if
                          (abs(x.location.start - dnstart) < 10)
                          and (abs(x.location.end - dnend) < 10)
                          and x.identifier.split(".")[0] == pfam.split(".")[0]]
                    if dn:
                        dn = dn[0]
                        ident_fn = lambda fident: "_".join(fident.split("_")[0:1] + fident.split("_")[-2:])
                        pdb = [x for x in protein.features if
                               x.type == feature_type and ident_fn(x.identifier) == ident_fn(hit.id)]
                        if not pdb:
                            posSet = set(range(dn.location.start, dn.location.end))
                            dncover = 1.0 * len(posSet & set(range(dnstart, dnend))) / (
                                    dn.location.end - dn.location.start)
                            if (min_identity <= identity(hsp)) and (dncover >= min_query_coverage):
                                hsp_feature = BioDocFactory.feature_from_hsp(hsp, feature_type)
                                hsp_feature.location.start += dnstart
                                hsp_feature.location.end += dnstart
                                features_added = features_added + 1
                                change = True
                                protein.features.append(hsp_feature)
                if change:
                    protein.save()

    _log.info("Features added: " + str(features_added))


def load_blast_features( organism, blast_file, feature_type,
                        min_identity=0, min_query_coverage=0, min_hit_coverage=0):
    queries = list(bpsio.parse(blast_file, 'blast-xml'))

    def check_overlap(features, new_feature, max_aa_overlap):
        for f in features:
            if (1.0 * len(new_feature & f) / len(f)) > 0.8:
                return True
        return False

    features_added = 0
    for query in tqdm(queries):

        gene = query.id

        proteins = Protein.objects(organism=organism, gene=gene).no_cache().timeout(False)
        change = False
        for protein in proteins:
            for hit in query:
                hsp = hit[0]
                if ((identity(hsp) >= min_identity)
                        and (coverage(query, hsp) >= min_query_coverage)
                        and (hit_coverage(hit, hsp) >= min_hit_coverage)):

                    hsp_feature = BioDocFactory.feature_from_hsp(hsp, feature_type)
                    features_added = features_added + 1
                    change = True
                    protein.features.append(hsp_feature)
                elif (identity(hsp) >= min_identity) and (hit_coverage(hit, hsp) >= min_hit_coverage):
                    for dn in protein.domains():
                        posSet = set(range(dn.location.start, dn.location.end))
                        dncover = 1.0 * len(posSet & set(range(hsp.query_start, hsp.query_end))) / (
                                dn.location.end - dn.location.start)
                        if dncover >= min_query_coverage:
                            hsp_feature = BioDocFactory.feature_from_hsp(hsp, feature_type)
                            features_added = features_added + 1
                            change = True
                            protein.features.append(hsp_feature)
            if change:
                protein.save()

    _log.info("Features added: " + str(features_added))
