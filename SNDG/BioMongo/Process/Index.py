"""

"""

import logging
import re

from tqdm import tqdm

from SNDG.BioMongo.Model.Ontology import Ontology
from SNDG.BioMongo.Model.Protein import Protein
from SNDG.BioMongo.Model.SeqColDruggabilityParam import SeqColDruggabilityParamTypes
from SNDG.BioMongo.Model.SeqColOntologyIndex import SeqColOntologyIndex
from SNDG.BioMongo.Model.SeqCollection import SeqCollection, Metric
from SNDG.BioMongo.Process.BioCyc2Mongo import BioCyc
from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
from SNDG.BioMongo.Process.EC2Mongo import EC2Mongo
from SNDG.BioMongo.Process.GO2Mongo import GO2Mongo
from SNDG.BioMongo.Process.KeywordIndexer import KeywordIndexer
from SNDG.BioMongo.Process.StructureIndexer import StructuromeIndexer
from functools import reduce

_log = logging.getLogger(__name__)

def index_seq_collection(db, genome, ec=True, go=True, keywords=True, organism_idx=True, pathways=True,structure=False):

    collection = SeqCollection.objects(name=genome).get()

    if ec:
        ec2mongo = EC2Mongo(db)
        _log.debug("Building EC index...")
        ec2mongo.pre_build_index(collection)
        collection.ec_index = True
        _log.debug("EC index finished")
        collection.save()

    if go:
        go2mongo = GO2Mongo("/data/databases/go/go.obo", db)
        go2mongo.init()
        _log.debug("Building GO index...")
        go2mongo.pre_build_index(collection)
        collection.go_index = True
        collection.save()
        _log.debug("GO index finished")



    if structure:
        si = StructuromeIndexer(collection)
        si.build_index()

    if pathways:
        biocyc = BioCyc(db)
        biocyc.user = BioMongoDB.demo
        _log.debug("Building Biocyc index...")
        biocyc.pre_build_index(collection)
        _log.debug("Biocyc index finished")

    if keywords:

        _log.debug("indexing by keyword...")
        ki = KeywordIndexer()
        cache = {}
        total_p = db.proteins.count({"organism":genome})
        with tqdm(Protein.objects(organism=genome).no_cache().timeout(False), total=total_p) as pbar:

            for prot in pbar:
                pbar.set_description( prot.name )
                # Basic keywords
                current_keywords = list(set([x.lower().strip() for x in reduce(list.__add__,
                                                                               map(ki.extract_keywords, [prot.name, prot.description] + prot.gene))]))

                prot.keywords = current_keywords + prot.keywords
                # ontologies keywords
                terms = prot.ontologies
                terms = terms + [x.identifier.strip().lower() for x in prot.features if x.identifier ]
                terms = terms + [x.type.strip().lower() for x in prot.features if x.identifier ]
                terms = list(set([x.lower() for x in terms]))

                for term in terms  :
                    if term not in cache:
                        ont = Ontology.objects(term=str(term))
                        if len(ont):
                            cache[term] = ont.first()

                    if term in cache:
                        prot.keywords = prot.keywords + cache[term].keywords
                    # SO:0001060 missense_variant

                prot.keywords = list(set(prot.keywords + terms))
                prot.save()
        _log.debug("Keyword index finished")



    if organism_idx:
        _log.debug("indexing ontology by organism")
        prots = list(db.proteins.find({"organism":genome, "ontologies.0":{"$exists":True}}))
        for prot in tqdm(prots):
            for term in prot["ontologies"]:
                if (term in cache) and cache[term].ontology not in ["ec", "go"]:
                    seq_col_ont_idx = SeqColOntologyIndex(seq_collection_id=collection.id, term=term,
                                                          seq_collection_name=genome, name=cache[term].name,
                                                          ontology=cache[term].ontology,
                                                          keywords=cache[term].keywords)
                    seq_col_ont_idx.save()
        SeqColOntologyIndex.objects(count=0).delete()
        _log.debug("Organism index finished")
    collection.save()
    _log.info("indexing %s finished" + genome)

def build_statistics(db, genome_name):
    import numpy as np
    genome = SeqCollection.objects(name=genome_name).get()  # @UndefinedVariable
    genome.statistics = []
    total = db.proteins.count({"organism": genome_name})

    drug = [x["search"]["druggability"] for x in db.proteins.find(
        {"organism": genome_name,
         "search.druggability": {"$exists": 1}},
        {"search.druggability": 1})]
    if drug:
        count, _ = np.histogram(drug, bins=10)
        freqs = 1.0 * count / sum(count)
        genome.statistics.append(Metric(name='druggabilityDistribution', values=list(freqs)))

    for dp in genome.druggabilityParams:
        if dp.target == "pathway":
            total = len(genome.pathways)
            if dp.type == SeqColDruggabilityParamTypes.number:
                data = [x.properties[dp.name]
                        for x in genome.pathways if dp.name in x.properties]
                count, division = np.histogram(data, bins=10)
                labels = [str(x) for x in division]
                count = [x for x in count]
            else:  # values
                labels = dp.options

                count = [len([1 for x in genome.pathways
                              if (dp.name in x.properties)
                              and (x.properties[dp.name] == label)])
                         for label in labels]
        else:
            total = db.proteins.count({"organism": genome_name, "search." + dp.name: {"$exists": 0}})
            if dp.type == SeqColDruggabilityParamTypes.number:
                data = [x["search"][dp.name]
                        for x in db.proteins.find(
                        {"organism": genome_name, "search." + dp.name: {"$exists": 1}},
                        {"search." + dp.name: 1}) if
                        isinstance(x["search"][dp.name], float) or isinstance(x["search"][dp.name], int)]
                count, division = np.histogram(data, bins=10)
                labels2 = [str(round(x, 2)) for x in division]
                labels = [labels2[i] + "-" + labels2[i + 1] for i, _ in list(enumerate(labels2))[:-1]]
                count = [x for x in count]
            else:  # values
                labels = dp.options

                count = [db.proteins.count({"organism": genome_name, "search." + dp.name:
                    (label if label not in ["true", "false"] else (label == "true"))})
                         for label in labels]

        labels = ["unannotated"] + labels

        genome.statistics.append(Metric(name='dp_' + dp.name, values=[total] + count, labels=labels))

    # Ontologies
    genome.statistics.append(
        Metric(name="proteins", value=db.proteins.find({"seq_collection_id": genome.id}).count()))
    ecre = re.compile("^ec:")
    genome.statistics.append(
        Metric(name="ec", value=db.proteins.find({"seq_collection_id": genome.id, "ontologies": ecre}).count()))
    gore = re.compile("^go:")
    genome.statistics.append(
        Metric(name="go", value=db.proteins.find({"seq_collection_id": genome.id, "ontologies": gore}).count()))

    genome.statistics.append(Metric(name="SO:0001079", description="polypeptide_structural_motif",
                                    value=db.proteins.find(
                                        {"organism": genome_name, "features.type": "SO:0001079"}).count()))

    #         db_structs = pymongo.MongoClient().pdb
    #         struct_count = list(db_structs.structures.aggregate([
    #             {"$match":{"organism":genome_name}}, {"$group":{"_id":"$templates.aln_query.name"}},
    #             {"$group":{"_id":"", "count":{"$sum":1}}}]))
    #         struct_count = struct_count[0]["count"] if struct_count else 0

    genome.statistics.append(Metric(name="Models", description="Structures generated by Homology",
                                    value=db.proteins.find(
                                        {"organism": genome_name, "search.structure_type": "model"}).count()));

    genome.statistics.append(
        Metric(name="SO:0000417", description="polypeptide_domain_count", value=db.proteins.find(
            {"organism": genome_name, "features.type": "SO:0000417"}).count()))
    genome.statistics.append(Metric(name="SO:0000418", description="signal_peptide_count", value=db.proteins.find(
        {"organism": genome_name, "features.type": "SO:0000418"}).count()))
    genome.statistics.append(
        Metric(name="SO:0100009", description="lipo_signal_peptide_count", value=db.proteins.find(
            {"organism": genome_name, "features.type": "SO:0100009"}).count()))
    genome.statistics.append(Metric(name="SO:0001077", description="transmembrane_count", value=db.proteins.find(
        {"organism": genome_name, "features.type": "SO:0001077"}).count()))

    genome.save()

if __name__ == '__main__':
    import pymongo
    from SNDG import init_log
    from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
    #from SNDG.BioMongo.Process.Index import index_seq_collection
    init_log()
    BioMongoDB("tdr")
    index_seq_collection(pymongo.MongoClient().tdr,"GCF_001624625.1",
                         ec=True, go=True, keywords=True, organism_idx=True, pathways=False,
                         structure=False)
