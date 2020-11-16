'''
Created on May 26, 2017

@author: esosa

hgvs.parser is used to parse http://varnomen.hgvs.org/

'''

import logging

import vcf
import hgvs.parser

_log = logging.getLogger(__name__)


class SnpeffEffect():
    hgvsparser = hgvs.parser.Parser()

    def __str__(self):
        return "Snpeff(%s; %s; %s; %s;)" % (
        self.geneid, str(self.annotation), str(self.hgvs_c), str(self.hgvs_p) if self.hgvs_p else "")

    def __repr__(self):
        return self.__str__()

    def __init__(self, alt, effects, impact, gene, geneid, feature_type, feature_id,
                 transcript_biotype, rank_div_total, hgvs_c, hgvs_p, c_dna_pos,
                 cds_pos, aa_pos, dist_to_feature, errors,
                 aa_len):
        self.alt = alt
        self.annotation = effects
        self.impact = impact
        self.gene = gene
        self.geneid = geneid
        self.feature_type = feature_type
        self.feature_id = feature_id
        self.transcript_biotype = transcript_biotype
        self.rank_div_total = rank_div_total
        self.hgvs_c = hgvs_c
        self.hgvs_p = hgvs_p
        self.c_dna_pos = c_dna_pos
        self.cds_pos = cds_pos
        self.aa_pos = aa_pos
        self.dist_to_feature = dist_to_feature
        self.errors = errors
        self.aa_len = aa_len
        self.aa_ref = ""
        self.aa_alt = ""
        if self.hgvs_c:
            self.gene_pos = hgvs_c.pos.start.base
        if self.hgvs_p and self.hgvs_p.pos:

            self.aa_ref = self.hgvs_p.pos.start.aa
            try:
                if self.hgvs_p.edit.type == "del":
                    self.aa_alt = "del"
                elif self.hgvs_p.edit.type == "dup":
                    self.aa_alt = "dup"
                else:
                    self.aa_alt = self.hgvs_p.edit.alt
            except:
                _log.warn(self.hgvs_p.edit)

    @classmethod
    def read(cls, ann_str):

        aa_pos, aa_len = (None, None)
        (alt, annotation, impact, gene, geneid, feature_type, feature_id, transcript_biotype,
         rank_div_total, hgvs_c, hgvs_p, c_dna_pos, cds_pos, (aa_pos_aa_len), dist_to_feature, errors) = ann_str.split(
            "|")
        annotation = annotation.split("&")
        if aa_pos_aa_len:
            aa_pos, aa_len = aa_pos_aa_len.split("/")
            aa_pos = int(aa_pos)
        try:
            hgvs_c = cls.hgvsparser.parse_hgvs_variant("xx:" + hgvs_c).posedit
        except Exception as ex:
            _log.warn(ex)
            hgvs_c = ""
        if hgvs_p:
            try:
                hgvs_p = cls.hgvsparser.parse_hgvs_variant("xx:" + hgvs_p).posedit
            except Exception as ex:
                _log.warn(ex)
                hgvs_p = ""
        else:
            hgvs_p = None

        return SnpeffEffect(alt, annotation, impact, gene, geneid, feature_type, feature_id, transcript_biotype,
                            rank_div_total, hgvs_c, hgvs_p, c_dna_pos, cds_pos, aa_pos, dist_to_feature, errors, aa_len)


class VcfSnpeffIO():

    @classmethod
    def parse(cls, vcf_path):
        if hasattr(vcf_path, "read"):
            h = vcf_path
        else:
            h = open(vcf_path)

        try:
            variantes = vcf.VCFReader(h)
            for v in variantes:
                effects = [SnpeffEffect.read(x) for x in (v.INFO["ANN"] if "ANN" in v.INFO else [])]
                intergenic = [(i, x) for i, x in enumerate(effects) if "intragenic_variant" in x.annotation]
                if intergenic:
                    i, intergenic = intergenic[0]
                    if (("upstream_gene_variant" in effects[0].annotation)
                            or ("downstream_gene_variant" in effects[0].annotation)):
                        effects = [effects[i]] + effects[:i - 1] + effects[i:]
                yield (v, effects)
        finally:
            h.close()
