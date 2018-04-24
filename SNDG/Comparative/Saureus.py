'''
Created on Oct 5, 2017

@author: eze
'''

import logging

import Bio.SeqIO as bpio
import mongoengine
import pandas as pd
from bson.objectid import ObjectId

from SNDG import init_log
from SNDG.BioMongo.Model.Feature import Feature, Location
from SNDG.BioMongo.Model.Protein import Protein
from SNDG.BioMongo.Model.SeqColDruggabilityParam import SeqColDruggabilityParam, SeqColDruggabilityParamTypes
from SNDG.BioMongo.Model.SeqCollection import SeqCollection

_log = logging.getLogger(__name__)


class Saureus(object):
    """
    http://mbio.asm.org/content/suppl/2016/05/03/mBio.00444-16.DCSupplemental -> Tabla  S1
    Manually curated  10/10/2017    
    """

    drugs = ["Penicillin", "Cefoxitin", "Ciprofloxacin", "Moxifloxacin", "Gentamicin", "Amikacin", "Tobramycin",
             "Erythromycin",
             "Clindamycin", "Tetracycline", "Fusidic_acid", "Linezolid", "Mupirocin", "Rifampicin", "Trimethoprim",
             "Daptomycin", "Vancomycin"
             ]

    def __init__(self, db, organism, user="demo"):
        self.db = db
        self.organism = organism
        self.user = user

    @staticmethod
    def parse_change(change):

        return (change[0], change[-1])

    def _process_prot(self, prot, r, i):

        try:
            pos = int(r.Substitution.split(",")[0][1:-1]) - 1
            start = pos
            end = pos
        except:
            if r.Substitution == "Deletions":
                start = 0
                end = len(prot.seq) - 1
            else:
                _log.warn("error parsing subtitution position: %s -> %s" % (r["Core gene"], r["Substitution"]))
                return

        quals = {
            "drug": r.Antibiotic,
            "change": r.Substitution,
            "gene": r["Core gene"]
        }
        if r.Reference:
            quals["reference"] = r.Reference

        fvariant = Feature(_id=ObjectId(), location=Location(start=start, end=end), type="Aanensen2016",
                           identifier="Aanensen2016_ " + str(i),
                           qualifiers=quals)
        prot.features.append(fvariant)
        prot.save()
        if self.user == "demo":
            prot.search.resistance = True
            prot.search[r.Antibiotic.lower()] = True
            prot.save()
        else:
            self.db.proteins.update({"_id": prot.id},
                {"$set": {"search." + self.user + ".resistance": True,
                          "search." + self.user + "." + r.Antibiotic.lower(): True}})



    def update_genome_props(self):
        if self.user == "demo":
            user2 = ""
        else:
            user2 = self.user + "."
        search_params = [("resistance", "Associated with resistance", "variant-db",
                          SeqColDruggabilityParamTypes.value, ["true", "false"], "true", "equal", "avg")

                         ]
        search_params = search_params + [
            (x.lower(), "Associated with " + x + " resistance", "variant-db",
             SeqColDruggabilityParamTypes.value, ["true", "false"], "true", "equal", "avg")
            for x in Saureus.drugs
        ]

        SeqCollection.objects(name=self.organism).update(
            __raw__={"$pull": {"druggabilityParams": {"target": "variant-db", "uploader": self.user}}})
        collection = SeqCollection.objects(name=self.organism).get()
        for name, description, target, _type, options, defaultValue, defaultOperation, defaultGroupOperation in search_params:
            Protein.objects(organism=self.organism).update(__raw__={"$set": {"search." + user2 + name: False}})
            if not collection.has_druggability_param(name):
                dp = SeqColDruggabilityParam(name=name, description=description, target=target,
                                             type=_type, uploader=self.user)
                dp.options = options
                dp.defaultValue = defaultValue
                dp.defaultOperation = defaultOperation
                dp.defaultGroupOperation = defaultGroupOperation
                collection.druggabilityParams.append(dp)
        collection.save()

    def create_fasta(self, data_csv="/data/projects/Staphylococcus/xomeq/resistance.csv"):
        seqs = open("/data/projects/Staphylococcus/annotation/curated/not_core.fasta", "w")
        uniprot = Uniprot()

        df = self.create_df(data_csv)
        try:
            for _, r in df.iterrows():
                if r.Obsolete:
                    continue
                gene = r["Accessory gene"] if r["Accessory gene"] else r["Core gene"]
                gene = gene.strip()
                assert gene

                unip = r["Protein Accession"]
                if unip:
                    seq = uniprot.download_and_load_seqrecord(unip)
                    seq.id = unip
                    seq.description = gene
                    bpio.write(seq, seqs, "fasta")
        finally:
            seqs.close()

    def create_df(self, data_csv="/data/projects/23staphylo/raw/metadata/saureus_resist_snps.csv"):
        df = pd.read_csv(data_csv)
        df["Substitution"] = df["Substitution"].str.strip()
        df["Antibiotic"] = df["Antibiotic"].str.strip()
        df["Core gene"] = df["Core gene"].str.strip()
        #        df["Accessory gene"] = df["Accessory gene"].str.strip()
        df = df.fillna(value="")
        return df

    def process(self, data_csv="/data/projects/23staphylo/raw/metadata/saureus_resist_snps.csv",
                genome="SaureusN315", user="demo"):
        Protein.objects(organism=genome).update(__raw__={"$pull": {"features": {"type": "Aanensen2016"}}})

        df = self.create_df(data_csv)
        for i, r in df.iterrows():

            gene = r["Core gene"]
            gene = gene.strip()
            assert gene

            # regex = re.compile(gene + '.*')
            has_prot = False
            for prot in Protein.objects(organism=genome, gene__iexact=gene):
                has_prot = True
                # print len([x.gene[1] if len(x.gene) > 1 else x.gene[0] for x in prot])
                self._process_prot(prot, r, i)
            if not has_prot:
                if "RNA" not in r["Core gene"]:
                    if r["Core gene"]:
                        print "Core gene not found %s" % gene
                    else:
                        print "%s not found" % gene


if __name__ == '__main__':
    init_log()
    from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB

    mdb = BioMongoDB("tdr")
    sa = Saureus(mdb.db, organism="SaureusN315", user="claudia")
    #     sa.create_fasta()
    sa.update_genome_props()
    sa.process(genome="SaureusN315", user="claudia")
