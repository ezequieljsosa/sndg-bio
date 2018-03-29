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
from SNDG.BioMongo.Model.SeqColDruggabilityParam import SeqColDruggabilityParamTypes, SeqColDruggabilityParam
from SNDG.BioMongo.Model.SeqCollection import SeqCollection
from SNDG.WebServices.Uniprot import Uniprot

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

    @staticmethod
    def parse_change(change):

        return (change[0], change[-1])

    def _process_prot(self, prot, r, i):
        mut = None
        try:
            pos = int(r.Substitution[1:-1])
        except:
            _log.warn("error parsing subtitution position: %s -> %s" % (r["Core gene"], r["Substitution"]))
            return

        quals = {
            "drug": r.Antibiotic,
            "change": r.Substitution,
            "gene": r["Core gene"],
            "uniprot": r["Protein Accession"]
        }
        if r.Reference:
            quals["reference"] = r.Reference
        if mut:
            quals["mut"] = mut
        fvariant = Feature(_id=ObjectId(), location=Location(start=pos, end=pos), type="Aanensen2016",
                           identifier="Aanensen2016_ " + str(i),
                           qualifiers=quals)
        prot.features.append(fvariant)
        prot.search.resistance = True
        prot.search[r.Antibiotic.lower()] = True
        prot.save()

    def update_genome_props(self, organism):
        search_params = [("resistance", "Associated with resistance", "variant-db",
                          SeqColDruggabilityParamTypes.value, ["true", "false"], "true", "equal", "avg")

                         ]
        search_params = search_params + [
            (x.lower(), "Associated with " + x + " resistance", "variant-db",
             SeqColDruggabilityParamTypes.value, ["true", "false"], "true", "equal", "avg")
            for x in Saureus.drugs
        ]

        SeqCollection.objects(name=organism).update(__raw__={"$pull": {"druggabilityParams": {"target": "variant-db"}}})
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

    def create_fasta(self, data_csv="saureus_resist_genes.csv", dst_seqs="resist.fasta"):
        uniprot = Uniprot()

        df = self.create_df(data_csv)
        with open(dst_seqs,"w") as h:
            for _, r in df.iterrows():
                # if r.Obsolete:
                #     continue
                gene = r["Accessory gene"]
                gene = gene.strip()

                unip = r["Protein Accession"]
                if unip:

                    seq = uniprot.download_and_load_seqrecord(unip)
                    if seq:
                        seq.id = unip
                        seq.description = gene
                        bpio.write(seq, h, "fasta")
                    else:
                        print unip

    def create_df(self, data_csv="/data/projects/Staphylococcus/xomeq/resistance.csv"):
        df = pd.read_csv(data_csv)
        #df["Substitution"] = df["Substitution"].str.strip()
        df["Antibiotic"] = df["Antibiotic"].str.strip()
        # df["Core gene"] = df["Core gene"].str.strip()
        df["Accessory gene"] = df["Accessory gene"].str.strip()
        #df = df.fillna(value="")
        return df

    def process(self, data_csv="/data/projects/Staphylococcus/xomeq/resistance.csv", genome="SaureusN315"):
        Protein.objects(organism=genome).update(__raw__={"$pull": {"features": {"type": "Aanensen2016"}}})

        df = self.create_df(data_csv)
        for i, r in df.iterrows():
            if r.Obsolete:
                continue
            gene = r["Accessory gene"] if r["Accessory gene"] else r["Core gene"]
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
    init_log(rootloglevel=logging.INFO)
    logging.getLogger("requests").setLevel(logging.WARN)
    #mongoengine.connect("saureus")
    sa = Saureus()

    sa.create_fasta("/data/projects/23staphylo/raw/metadata/saureus_resist_genes.csv",
                    "/data/projects/23staphylo/processed/resist/resistdb.fasta")
    #sa.update_genome_props(organism="SAureusN315")
    #sa.process(genome="SAureusN315")
    print ("OK")
