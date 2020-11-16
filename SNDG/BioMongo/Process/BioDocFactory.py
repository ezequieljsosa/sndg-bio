'''
Created on May 22, 2017

@author: eze
'''

import zlib
from bson.objectid import ObjectId
from Bio.Seq import Seq
from SNDG.BioMongo.Model.SeqCollection import TaxEDoc
from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
from SNDG.BioMongo.Model.Sequence import BioProperty, Size, Contig
from SNDG.BioMongo.Model.Feature import Feature, Location
from SNDG.BioMongo.Model.Protein import Protein
from SNDG.BioMongo.Model.SeqCollection import Genome

from SNDG.BioMongo.Model.Alignment import SimpleAlignment, AlnLine
from SNDG.WebServices.NCBI import NCBI
from SNDG import Struct


class BioDocFactory(object):
    '''
    BioDocFactory
    '''

    @staticmethod
    def create_genome(name, seqrecord, ncbi_tax=None, tax_provider=None):
        '''
        TaxEDoc
        '''
        tax = None
        if tax_provider:
            if not ncbi_tax:
                ncbi_tax = seqrecord.annotations["ncbi_taxid"] if "ncbi_taxid" in seqrecord.annotations else None

            tax = tax_provider.getTax(ncbi_tax) if ncbi_tax else None
            tax = TaxEDoc(tid=tax.ncbi_taxon_id, superkingdom=tax.node_rank,
                          name=[x.name for x in tax.names if x.name_class == tax_provider.scientific_name][
                              0]) if tax else None
        return Genome(
            name=name,
            description=seqrecord.description,
            organism=seqrecord.annotations["organism"] if "organism" in seqrecord.annotations else
            seqrecord.annotations["taxonomy"][-1] if "taxonomy" in seqrecord.annotations else "",
            auth=str(BioMongoDB.demo_id),
            tax=tax
        )

    @staticmethod
    def create_features_from_contig(seqrecord, source, type_map={x: x for x in NCBI.ftypes},
                                    extract_annotation_feature=lambda feature: feature):
        ftypes = {xx: 1 for xx in type_map}
        features = []
        gene_ids = {}
        for feature in seqrecord.features:
            f = extract_annotation_feature(feature)
            if f.type in ftypes:

                fid = f.qualifiers["description"][0] if "description" in f.qualifiers else f.id
                if "ID" in f.qualifiers:
                    fid = f.qualifiers["ID"][0]
                if "product" in f.qualifiers:
                    fid = f.qualifiers["product"][0]
                if "gene_id" in f.qualifiers:
                    fid = f.qualifiers["gene_id"][0]
                if "gene" in f.qualifiers:
                    fid = f.qualifiers["gene"][0]
                if "protein_id" in f.qualifiers:
                    fid = f.qualifiers["protein_id"][0]
                if "tRNA_anti-codon" in f.qualifiers:
                    fid = fid + " -> " + f.qualifiers["tRNA_anti-codon"][0]

                fdoc = Feature(_id=ObjectId(), identifier=fid,
                               location=Location(start=f.location.start,
                                                 end=f.location.end,
                                                 strand=f.location.strand),
                               type=type_map[feature.type])

                if "locus_tag" in f.qualifiers:
                    locus_tag = f.qualifiers["locus_tag"][0]
                    fdoc.identifier = locus_tag
                    fdoc.locus_tag = locus_tag
                    fdoc.alias.append(fdoc.locus_tag)
                else:
                    fdoc.locus_tag = fid

                gene_ids[fdoc.locus_tag] = fdoc._id
                if "gene" in f.qualifiers:
                    fdoc.alias.append(f.qualifiers["gene"][0])
                if "protein_id" in f.qualifiers:
                    fdoc.alias.append(f.qualifiers["protein_id"][0])

                if "old_locus_tag" in f.qualifiers:
                    fdoc.alias = fdoc.alias + f.qualifiers["old_locus_tag"]
                if source:
                    fdoc.source = source
                features.append(fdoc)
        return (features, gene_ids)

    @staticmethod
    def create_contig(seqrecord, seq_col, source=None, type_map={x: x for x in NCBI.ftypes},
                      extract_annotation_feature=lambda feature: feature
                      ):

        features, gene_ids = BioDocFactory.create_features_from_contig(
            seqrecord, source, type_map, extract_annotation_feature=extract_annotation_feature,
        )

        seq = str(seqrecord.seq)
        bigseq = None
        if len(seq) > ((1024 ** 2) * 12):
            bigseq = zlib.compress(str(seqrecord.seq).encode('utf-8'))

        contig = Contig(
            organism=seq_col.name,
            seq_collection_id=seq_col.id,
            name=seqrecord.id,
            description=seqrecord.description,
            seq="",
            features=features,
            properties=[BioProperty(_type="annotation",
                                    property=xx, value=seqrecord.annotations[xx][0])
                        for xx in ["topology", "molecule_type"]
                        if xx in seqrecord.annotations],

            size=(Size(unit="bp", len=len(str(seqrecord.seq)))),
            auth=str(BioMongoDB.demo_id))
        if bigseq:
            contig.bigseq = bigseq
        else:
            contig.seq = seq

        return contig, gene_ids

    @staticmethod
    def feature_from_hsp(hsp, feature_type):
        return Feature(_id=ObjectId(), location=Location(start=hsp.query_start, end=hsp.query_end),
                       aln=SimpleAlignment(
                           evalue=hsp.evalue,
                           aln_query=AlnLine(name=hsp.query.id, seq=str(hsp.aln[0].seq), start=hsp.query_start,
                                             end=hsp.query_end),
                           aln_hit=AlnLine(name=hsp.hit.id, seq=str(hsp.aln[1].seq), start=hsp.hit_start,
                                           end=hsp.hit_end),
                           aln_mid=hsp.aln_annotation["similarity"] if "similarity" in hsp.aln_annotation else ""
                       ),
                       identifier=hsp.hit.id, type=feature_type)

    @classmethod
    def alias(cls, feature):
        alias = [feature.id]

        if "gene_symbol" in feature.qualifiers:
            alias += feature.qualifiers["gene_symbol"]
        if "old_locus_tag" in feature.qualifiers:
            alias += feature.qualifiers["old_locus_tag"]
        if "locus_tag" in feature.qualifiers:
            alias += feature.qualifiers["locus_tag"]
        if "protein_id" in feature.qualifiers:
            alias += feature.qualifiers["protein_id"]
        if "Alias" in feature.qualifiers:
            alias += feature.qualifiers["Alias"]
        return [x for x in (set([x for x in set(alias) if x and ("unknown" not in x)])) if x.strip()]

    @classmethod
    def create_protein(cls, seqrecord, feature, exons=[]):

        alias = cls.alias(feature)

        locus_tag = feature.qualifiers["locus_tag"][0]
        protein_name = locus_tag
        seq = str(seqrecord.seq)
        # if "translation" in feature.qualifiers:
        #     seq = feature.qualifiers["translation"][0]
        # else:
        #     if exons:
        #         seq = str(reduce(Seq.__add__, [exon.extract(seqrecord.seq) for exon in exons]).translate())
        #     else:
        #         seq = str(feature.extract(Seq(str(seqrecord.seq))).translate())

        if "gene_symbol_source" in feature.qualifiers:
            try:
                int(feature.qualifiers["gene_symbol_source"][0])
                protein_name = ""
            except:
                protein_name = feature.qualifiers['gene_symbol'][0]
        #         elif "product" in mrna_feature.qualifiers:
        #             protein_name = mrna_feature.qualifiers['product'][0]

        p = Protein(seq=seq, name=protein_name)

        if "description" in feature.qualifiers:
            protein_description = feature.qualifiers['description'][0]
        elif "Note" in feature.qualifiers:
            protein_description = feature.qualifiers['Note'][0]
        elif "product" in feature.qualifiers:
            protein_description = feature.qualifiers['product'][0]
        else:
            protein_description = ""

        bp = BioProperty(_type="annotation", description="homolog proteins and sources used for the annotation")

        if "protein_id" in feature.qualifiers:
            bp.ncbi_protein_id = feature.qualifiers["protein_id"][0]

        if "db_xref" in feature.qualifiers:
            bp.ncbi_db_xref = feature.qualifiers["db_xref"][0]

        if "Dbxref" in feature.qualifiers:
            bp.ncbi_db_xref = feature.qualifiers["Dbxref"][0]

        if "top_cog_hit" in feature.qualifiers:
            bp.cog = feature.qualifiers["top_cog_hit"][0]
        if "gene_symbol_source" in feature.qualifiers:
            bp.source = feature.qualifiers["gene_symbol_source"][0]
        if "gene_product_name_source" in feature.qualifiers:
            bp.source = feature.qualifiers["gene_product_name_source"][0]
        ecs = []

        if "EC" in feature.qualifiers:
            ecs = ecs + [x.lower() if x.lower().startswith("ec") else "ec:" + x.lower() for x in
                         feature.qualifiers["EC"] if "." in x]

        if "EC_number" in feature.qualifiers:
            ecs = ecs + [x.lower() if x.lower().startswith("ec") else "ec:" + x.lower() for x in
                         feature.qualifiers["EC_number"] if "." in x]

        if "Dbxref" in feature.qualifiers:
            ecs = ecs + [x.lower() if x.lower().startswith("ec") else "ec:" + x.lower() for x in
                         feature.qualifiers["Dbxref"] if "." in x]

        gos = []

        if "db_xref" in feature.qualifiers:
            gos = gos + [x.lower() for x in feature.qualifiers["db_xref"] if
                         "GO:" in x and (x not in ["GO:0008150", "GO:0003674", "GO:0005575"])]
        if "GO" in feature.qualifiers:
            gos = gos + [x.lower() for x in feature.qualifiers["GO"] if
                         "GO:" in x and (x not in ["GO:0008150", "GO:0003674", "GO:0005575"])]
        if "GO_terms" in feature.qualifiers:
            gos = gos + [x.lower() for x in feature.qualifiers["GO_terms"] if
                     "GO:" in x and (x not in ["GO:0008150", "GO:0003674", "GO:0005575"])]
        if "Ontology_term" in feature.qualifiers:
            gos = [x.lower() for x in feature.qualifiers["Ontology_term"] if
                   "GO:" in x and (x not in ["GO:0008150", "GO:0003674", "GO:0005575"])]

        ontologies = list(set(ecs + gos))

        p.gene = list([locus_tag, protein_name]) if locus_tag != protein_name else [locus_tag]

        p.name = protein_name
        p.description = protein_description
        p.ontologies = ontologies
        p.properties = [bp]

        p.alias = alias

        return p
