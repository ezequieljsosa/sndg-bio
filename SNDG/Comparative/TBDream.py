'''
Created on Aug 11, 2014

@author: eze
'''
import logging
import math
from collections import defaultdict

import Bio.SeqUtils as bpsutils
import pandas as pd
from Bio import SeqUtils
from tqdm import tqdm
from SNDG.Comparative.VcfSnpeffIO import VcfSnpeffIO

_log = logging.getLogger(__name__)


class TBDream:
    '''
    Base de SNPs de Resistencia
    https://tbdreamdb.ki.se/Info/
    '''

    """
    rna positions in NC_000962.3
    """
    rna_pos = {'RVnc0001': 2139418,
               'RVnc0002': 918263,
               'RVnc0003': 2299744,
               'RVnc0004': 4099385,
               'RVnc0005': 704186,
               'RVnc0006': 4168153,
               'RVnc0008': 293603,
               'RVnc0010': 1914961,
               'RVnc0012': 1283692,
               'RVnc0013': 1413093,
               'RVnc0015': 1535416,
               'RVnc0018': 2517031,
               'RVnc0021': 1471618,
               'RVnc0024': 2692171,
               'RVnc0034': 1220387,
               'RVnc0035': 1547128,
               'RVnc0036': 1960666,
               'RVnc0036a': 4100668,
               'RVnc0040': 4317072,
               'RVnc0046': 3467966,
               'RVnc0047': 1952290,
               'Rvnr01': 1471845,
               'Rvnr02': 1473657,
               'Rvnr03': 1476898,
               'Rvns01': 2500444,
               'Rvnt01': 10886,
               'Rvnt02': 11111,
               'Rvnt03': 25643,
               'Rvnt04': 386203,
               'Rvnt05': 658108,
               'Rvnt06': 731493,
               'Rvnt07': 731602,
               'Rvnt08': 733523,
               'Rvnt09': 850641,
               'Rvnt10': 923802,
               'Rvnt11': 923998,
               'Rvnt12': 924109,
               'Rvnt13': 924212,
               'Rvnt14': 1025320,
               'Rvnt15': 1113510,
               'Rvnt16': 1138075,
               'Rvnt17': 1177395,
               'Rvnt18': 1446192,
               'Rvnt19': 1512727,
               'Rvnt20': 1828014,
               'Rvnt21': 1946612,
               'Rvnt22': 2401986,
               'Rvnt23': 2510597,
               'Rvnt24': 2581763,
               'Rvnt25': 2619406,
               'Rvnt26': 2765330,
               'Rvnt27': 2765540,
               'Rvnt28': 2794175,
               'Rvnt29': 2827853,
               'Rvnt30': 2835493,
               'Rvnt31': 2969496,
               'Rvnt32': 2969752,
               'Rvnt33': 2969854,
               'Rvnt34': 2969941,
               'Rvnt35': 3348546,
               'Rvnt36': 3348658,
               'Rvnt37': 3431839,
               'Rvnt38': 3559369,
               'Rvnt39': 4081364,
               'Rvnt40': 4126540,
               'Rvnt41': 4168344,
               'Rvnt42': 4199130,
               'Rvnt43': 4216864,
               'Rvnt44': 4216967,
               'Rvnt45': 4222580
               }

    drugs = ["AMI", "EMB", "FLQ", "INH", "PAS", "PZA", "RIF", "SM", "OTH"]

    def __init__(self):
        self.csv_db_path = "/data/databases/tbdream/" + "DownloadDB_CSV.csv"
        self._df = None

    def _codon_position(self, dict_snp):
        try:
            return int(dict_snp["EstimatedCodonPosition"])
        except:
            try:
                return int(dict_snp["ReportedCodonPosition"])
            except:
                try:
                    return (int(dict_snp["NucleotidePosition"]) / 3) + 1
                except:
                    # Ej: 1213-1214
                    try:
                        return (int(dict_snp["NucleotidePosition"].split("-")) / 3) + 1
                    except:
                        return None

    def _mutated_aa(self, dict_snp):
        '''
        returns tuple (aa_original,aa_replacement) or none if it can't be parsed    
        '''
        aminoacid = dict_snp["AminoAcid"].replace(" /", "/").replace("/ ", "/")
        aa = aminoacid.split(" ")[0]
        if "/" in aa:
            (aa_original, aa_replacement) = aa.split("/")
            if aa_original != "STOP" and aa_replacement != "STOP" and len(aa_replacement) == 3 and len(
                    aa_original) == 3:
                return (aa_original, aa_replacement)
            else:
                return (aa_original, aa_replacement)
        return None

    def _mutation_type(self, dict_snp):
        '''
        returns tuple (aa_original,aa_replacement) or none if it can't be parsed    
        '''
        aminoacid = dict_snp["AminoAcid"].replace(" /", "/").replace("/ ", "/")
        aa = aminoacid.split(" ")[0]
        if "Frameshift" in aa:
            return "frameshift"
        if "/" in aa:
            (aa_original, aa_replacement) = aa.split("/")
            if aa_original != "STOP" and aa_replacement != "STOP" and len(aa_replacement) == 3 and len(
                    aa_original) == 3:
                return "missense"
            elif aa_original == "STOP":
                return "stoploss"
            elif aa_replacement == "STOP":
                return "stopgain"
            else:
                return "unknown"
        return None

    def nu_pos(self, nu_pos):
        try:
            return int(nu_pos.split("-")[0])
        except:
            return -1

    def nu_alt(self, x):
        try:
            x.split("/")[1].strip()
        except:
            return "?"

    def load(self):
        '''
            GeneID    Drug    GeneName    NucleotidePosition    Polymorphism    EstimatedCodonPosition    ReportedCodonPosition    AminoAcid
        0    Rv3795    EMB    embB        1121-1122            GGC/GTG            374                        374                    Gly/Val
        1    Rv3795    EMB    embB        1123-1124            CCG/GCG            375                        375                    Pro/Ala
        '''
        cols = ['ID', 'GeneID', 'SeqNo', 'Drug', 'GeneName', 'ApprovalLevel', 'AddedByUser', 'PrimaryReference',
                'ReviewReferece',
                'SecondaryReference', 'NucleotidePosition', 'Polymorphism', 'EstimatedCodonPosition',
                'ReportedCodonPosition',
                'AminoAcid', 'Note', 'TimePeriod', 'StudyPopulation', 'Country', 'MolecularDetectionMethod',
                'GeneCoverage',
                'ResistancePattern', 'MIC', 'SusceptibilityTestingMethod', 'RTotalIsolates', 'RSIsolatesWMutation',
                'AdditionalMutations',
                'HighQuality', 'PMID', 'OverallSequentialNumber', 'NoHC', 'Temp1', 'Temp2', 'Temp3', 'Temp4', 'Temp5',
                'MutationType']
        self._df = pd.DataFrame()

        errors = []
        i = 0
        rvs = []
        with open(self.csv_db_path, "r") as handle:
            lines = handle.readlines()
            _log.debug("Total lines count: " + str(len(lines)))
            for idx, line in enumerate(lines):
                i = i + 1
                if i < 3:
                    continue

                split = line.split(",")
                if len(split) >= 14:

                    dict_snp = {"line": line}
                    j = 0
                    for col in cols:
                        dict_snp[col] = split[j]
                        j = j + 1
                    if (not (("/" in dict_snp["Polymorphism"]) or ("ins" in dict_snp["Polymorphism"]) or (
                            "del" in dict_snp["Polymorphism"]))
                            and
                            not (("/" in dict_snp["AminoAcid"]) or ("ins" in dict_snp["AminoAcid"]) or (
                                    "del" in dict_snp["AminoAcid"]))
                    ):
                        print dict_snp["ID"]

                    self._df = self._df.append(dict_snp, ignore_index=True)
                    rvs.append(dict_snp["GeneID"])
                    # i = i -1
                    # if not i:
                    #    break
                    # (rv,nucleotide,polimorphism,codon_aa,change_aa,drug) =  (split[1],split[10],split[11],split[14],split[12],split[3])
                    # if isinstance( codon_aa,int):
                    #    snps.append((rv,nucleotide,polimorphism,codon_aa,change_aa,drug))
                    # else:
                    #    errors.append(line)
                else:
                    errors.append(line)

        self._df["codon"] = map(self._codon_position, [x for i, x in self._df.iterrows()])
        self._df["change"] = map(self._mutated_aa, [x for i, x in self._df.iterrows()])

        self._df["ref"] = map(lambda x: bpsutils.seq1(x[0]) if x else None, self._df["change"])
        self._df["mutation"] = map(lambda x: bpsutils.seq1(x[1]) if x else None, self._df["change"])
        self._df["mut_type"] = map(self._mutation_type, [x for i, x in self._df.iterrows()])

        self._df["rv"] = map(lambda x: x.lower(), self._df["GeneID"])
        self._df["gene"] = map(lambda x: x.lower(), self._df["GeneName"])
        self._df["nu_pos"] = map(self.nu_pos, self._df["NucleotidePosition"])
        self._df["nu_ref"] = map(lambda x: x.split("/")[0].strip(), self._df["Polymorphism"])
        self._df["nu_alt"] = map(self.nu_alt, self._df["Polymorphism"])
        self._df["raw"] = line
        self._df["rna"] = map(lambda x: True if x == None else False, self._df["change"])

        _log.info("SNPs loaded:" + str(len(self._df)))
        _log.info("Errors: " + str(len(errors)))
        _log.info("RV count: " + str(len(set(rvs))))

    def __iter__(self):
        '''
        returns iterator of i,row
        row fields: rv, gene, ref, mutation, change
            and 'ID', 'GeneID', 'SeqNo', 'Drug', 'GeneName', 'ApprovalLevel', 'AddedByUser', 'PrimaryReference', 'ReviewReferece', 'SecondaryReference', 'NucleotidePosition', 'Polymorphism', 'EstimatedCodonPosition', 'ReportedCodonPosition', 'AminoAcid', 'Note', 'TimePeriod', 'StudyPopulation', 'Country', 'MolecularDetectionMethod', 'GeneCoverage', 'ResistancePattern', 'MIC', 'SusceptibilityTestingMethod', 'RTotalIsolates', 'RSIsolatesWMutation', 'AdditionalMutations', 'HighQuality', 'PMID', 'OverallSequentialNumber', 'NoHC', 'Temp1', 'Temp2', 'Temp3', 'Temp4', 'Temp5', 'MutationType']
        '''
        return self._df.iterrows()

    def process_vcf(self, vcf):
        data = defaultdict(lambda: defaultdict(lambda: {}))
        for variant, effects in tqdm(VcfSnpeffIO.parse(vcf)):
            effect = effects[0]
            for sample in variant.samples:
                sample_name = sample.sample.split(".variant")[0]
                if sample.called and "synonymous_variant" not in effect.annotation:

                    if self.variant_rv(effect.geneid):
                        level = "gene"
                        drugs = set([x["drug"] for x in self.variant_rv(effect.geneid)])
                        for drug in drugs:
                            data[sample_name][effect.geneid][drug] = level
                        pos = effect.aa_pos if effect.geneid not in self.rna_pos else (
                                effect.gene_pos - self.rna_pos[effect.geneid])
                        if (pos != None) and self.variant_pos(effect.geneid, pos):
                            level = "pos"
                            drugs = set([x["drug"] for x in self.variant_pos(effect.geneid, pos)])
                            for drug in drugs:
                                data[sample_name][effect.geneid][drug] = level

                            alt = str(variant.ALT[int(sample.data.GT) - 1])
                            ref = effect.aa_ref if effect.geneid not in self.rna_pos else variant.REF
                            alt = effect.aa_alt if effect.geneid not in self.rna_pos else alt
                            if self.exists_variant(effect.geneid, pos, ref, alt):
                                level = "reported"
                                drug = self.exists_variant(effect.geneid, pos, ref, alt)
                                data[sample_name][effect.geneid][drug] = level

        return {x: dict(y) for x, y in data.items()}

    def variant_rv(self, rv):
        row = self._df[(self._df["rv"] == rv.lower()) | (self._df["gene"] == rv.lower())]
        if len(row):
            return [{"ref": x["ref"], "change": x["mutation"], "pos": x["codon"], "drug": x["Drug"],
                     "pattern": x["ResistancePattern"]} for i, x in row.iterrows()]
        return None

    def variant_pos(self, rv, aa_pos):
        '''
        Supone aa_pos numerado desde 1
        '''
        row = self._df[(self._df["rv"] == rv.lower()) & ((self._df["codon"] == int(aa_pos)) |
                                                         (self._df["rna"] & (self._df["nu_pos"] == int(aa_pos))))]
        if len(row):
            return [{"ref": x["ref"], "change": x["mutation"], "pos": x["codon"], "drug": x["Drug"],
                     "pattern": x["ResistancePattern"]} for i, x in row.iterrows()]
        return None

    def variant_pos_df(self, rv, aa_pos):
        '''
        Supone aa_pos numerado desde 1
        '''
        row = self._df[(self._df["rv"] == rv.lower()) & (self._df["codon"] == int(aa_pos))]
        if len(row):
            return row.to_dict()
        return {}

    def exists_variant(self, rv, aa_pos, aa_ref, aa_change):
        '''
        Supone aa_pos numerado desde 1
        '''
        row = self._df[(self._df["rv"] == rv.lower()) & (self._df["codon"] == int(aa_pos))
                       & (self._df["ref"] == aa_ref) & (self._df["mutation"] == aa_change)]
        if len(row):
            return reduce(lambda x, y: " ".join([x, y]), set([x["Drug"] for i, x in row.iterrows()]))
        return None

    def parse_change(self, change):
        dref, dalt = change.split("->")[0].strip().split("/")
        dref = SeqUtils.seq1(dref)
        dalt = SeqUtils.seq1(dalt)
        return (dref, dalt)

    def load_in_sndg(self, organism="H37Rv"):
        from SNDG.BioMongo.Model.Protein import Protein
        from SNDG.BioMongo.Model.Feature import Feature, Location
        from SNDG.BioMongo.Model.SeqCollection import SeqCollection
        from SNDG.BioMongo.Model.SeqColDruggabilityParam import SeqColDruggabilityParamTypes, SeqColDruggabilityParam

        from bson.objectid import ObjectId

        search_params = [("resistance", "Associated with resistance", "variant-db",
                          SeqColDruggabilityParamTypes.value, ["true", "false"], "true", "equal", "avg")

                         ]
        search_params = search_params + [
            (x, "Associated with " + x + " resistance", "variant-db",
             SeqColDruggabilityParamTypes.value, ["true", "false"], "true", "equal", "avg")
            for x in TBDream.drugs
        ]

        Protein.objects(organism=organism).update(__raw__={"$pull": {"features": {"type": "tbdream"}}})
        collection = SeqCollection.objects(name=organism).get()
        for name, description, target, _type, options, defaultValue, defaultOperation, defaultGroupOperation in search_params:
            Protein.objects(organism=organism).update(__raw__={"$set": {"search." + name: False}})
            if not collection.has_druggability_param(name):
                dp = SeqColDruggabilityParam(name=name, description=description, target=target,
                                             type=_type, uploader="demo")
                dp.options = options
                dp.defaultValue = defaultValue
                dp.defaultOperation = defaultOperation
                dp.defaultGroupOperation = defaultGroupOperation
                collection.druggabilityParams.append(dp)
        collection.save()

        for rv, rows in self._df.groupby("rv"):
            prot = list(Protein.objects(organism=organism, gene__iexact=rv))
            if prot:
                prot = prot[0]
                for _, r in rows.iterrows():
                    mut = None
                    if r.change:
                        change = str(r.change[0]) + "/" + str(r.change[1])
                        mut = SeqUtils.seq1(r.change[1])
                    else:
                        change = r.AminoAcid
                    if math.isnan(r.codon):
                        try:
                            pos = int(r.AminoAcid)
                        except:
                            _log.warn("couldnt find the variant position")
                            continue
                    else:
                        pos = int(r.codon)

                    try:
                        res, t = r.RTotalIsolates.strip().split("/")
                        r_div_total_coef = int(res) * 1.0 / int(t)
                        r_div_total = r.RTotalIsolates.strip()

                    except:
                        r_div_total = None
                        r_div_total_coef = None

                    quals = {
                        "drug": r.Drug,
                        "change": change,
                        "gene": r.GeneID,
                        "pattern": r.ResistancePattern,
                        "additional": r.AdditionalMutations,
                        "r_div_total": r_div_total,
                        "r_div_total_coef": r_div_total_coef,
                        "mic": r.MIC}
                    if mut:
                        quals["mut"] = mut
                    fvariant = Feature(_id=ObjectId(), location=Location(start=pos, end=pos), type="tbdream",
                                       identifier="TBDream id " + r.ID,
                                       qualifiers=quals)
                    prot.features.append(fvariant)
                    prot.search.resistance = True
                    prot.search[r.Drug] = True
                prot.save()


if __name__ == '__main__':
    from SNDG.Comparative.TBDream import TBDream

    tbd = TBDream()
    tbd.csv_db_path = "/home/eze/Downloads/DownloadDB_CSV_corregido_28feb2018.csv"
    tbd.load()
    tbd.process_vcf("/data/projects/mtbxdr/processed/h37rv_processing/strains.ann.gvcf")
