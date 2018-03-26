'''
Created on Sep 16, 2016

@author: eze
'''

import logging
import pymongo
from mongoengine.connection import connect
from pymongo.mongo_client import MongoClient

from SNDG import init_log
from SNDG.BioMongo.Model.SeqCollection import Genome, SeqColDruggabilityParam
from SNDG.BioMongo.Process.BioCyc2Mongo import BioCyc
from SNDG.BioMongo.Process.StructureAnotator import StructureAnotator
from SNDG.BioMongo.Process.StructureIndexer import StructuromeIndexer
from SNDG.BioMongo.Process.PathwaysAnnotator import PathwaysAnnotator
from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
from pprint import pprint

init_log("/tmp/validate.log")
_log = logging.getLogger(__name__)
db = MongoClient().tdr

import re

regx_go = re.compile("^go:")  # , re.IGNORECASE)
regx_ec = re.compile("^ec:")


def validate_pathways_protein(g):
    for dp in BioCyc.pathways_search_params:
        if not g.has_druggability_param(dp[0]):
            print "%s no tiene el atributo %s" % (g.name, dp[0])

    count = db.proteins.count({"organism": g.name, "reactions.0": {"$exists": 1}})
    if count == 0:
        print  g.name + " NO TIENE PROTEINAS CON REACCIONES!!! "

    centralidades = [x["search"]["centrality"] for x in db.proteins.find({"organism": g.name}, {"search": 1}) if
                     "centrality" in x["search"]]
    max_c = max(centralidades)
    if max_c != 1:
        print "hay un error en la centralidad cargada " + str(max_c) + " " + g.name

    chokepoints = db.proteins.count({"organism": g.name, "search.chokepoint": True})

    if chokepoints < 10:
        print " hay muy pocos chokepoints en " + g.name


def validate_pathways(g):
    # PathwaySumary
    #         term = StringField( required=True)
    #     name = StringField( required=True)
    #     count = IntField( required=True)

    dps = {dp[0]: [] for dp in BioCyc.pathways_search_params}

    for ps in g.pathways:
        for dp in dps:
            dps[dp].append(ps.properties[dp])

    for dp, values in dps.items():

        try:
            suma = sum(values)
            if dp == "max_centrality":
                centralities = values
            if suma == 0:
                print "no hay " + dp + " en los pws de " + g.name
        except:
            if len(set(values)) == 1:
                print "solo hay un valor de " + dp + "(" + values[0] + ") en los pws de " + g.name

    if (max(centralities) < 0.99):
        print "mal la centralidad en los pws de " + g.name + " : " + str(max(centralities))

    validate_pathways_protein(g)


def validate_prots(g):
    gos = db.proteins.count({"organism": g.name, "ontologies": regx_go})
    if gos < 1000:
        print "%s tiene pocas proteinas con go anotados: %i" % (g.name, gos)
    ecs = db.proteins.count({"organism": g.name, "ontologies": regx_ec})
    if ecs < 1000:
        print "%s tiene pocas proteinas con ec anotados: %i" % (g.name, ecs)


def check_prot_params(g, prot_params):
    for dp in prot_params:
        if "overexpression" in dp[0]:
            continue
        if "essentiality" in dp[0]:
            continue
        if not g.has_druggability_param(dp[0]):
            print "%s no tiene el atributo %s" % (g.name, dp[0])
        else:
            if dp[3] == "number":
                param_count = db.proteins.count({"organism": g.name, "search." + dp[0]: {"$gt": 0}})  # {"$exists":True}

            elif dp[4] and "true" not in dp[4]:
                param_count = db.proteins.count({"organism": g.name, "search." + dp[0]: {"$exists": True}})
            else:
                param_count = db.proteins.count({"organism": g.name, "search." + dp[0]: True})
            if not param_count:
                print "%s no tiene el proteinas con el atributo %s" % (g.name, dp[0])


def validate_protein_search(g):
    check_prot_params(g, BioCyc.protein_pathway_search_params)
    # check_prot_params(g, SeqColDruggabilityParam.default_params)
    # check_prot_params(g, StructuromeIndexer.search_params)


def validate_genome(g):
    validate_pathways(g)
    validate_protein_search(g)

    for x in ["ec", "go", "biocyc_pw", "biocyc_reac"]:
        if db.col_ont_idx.count({"ontology": x, "seq_collection_name": g.name}) == 0:
            print  g.name + " sin indice de organismos " + x


#     clean no chokes
#     db.proteins.update({"search.chokepoint":False,"properties.property" : "chokepoint"},
#                                         {  "$pull":{"properties":{"property":"chokepoint"}} },multi=True
#                                         ) 

#     clean no structure
#     sdps = [x[0] for x in StructuromeIndexer.search_params   if (x[0] != "has_structure")]
#     q = {"organism":g.name,"search.has_structure":False, 
#          "$or":[{ "search." + x: {"$exists":True} } for x in sdps  ]
#          
#     }
#     print db.proteins.count(q    )
#     unset = {k:"" for k in sdps}
#     unset["search.structures"] = "" 
#     print db.proteins.update( q, {"$unset":unset},multi=True )


if __name__ == '__main__':
    connect("tdr")
    genomes = list(Genome.objects(auth=BioMongoDB.demo_id))

    assert 13 == len(genomes), len(genomes)

    for genome in genomes:
        dps = [dp[0] for dp in
               SeqColDruggabilityParam.default_params + StructuromeIndexer.search_params +
               BioCyc.protein_pathway_search_params + BioCyc.pathways_search_params]
        genome.druggabilityParams = [x for x in genome.druggabilityParams if x.name not in dps]
        for name, description, target, _type, options, _, _, _ in (
                BioCyc.protein_pathway_search_params + BioCyc.pathways_search_params):
            dp = SeqColDruggabilityParam(name=name, description=description, target=target,
                                         type=_type, uploader="demo")
            genome.druggabilityParams.append(dp)

    biocyc = BioCyc(db)
    biocyc.user = "demo"

    mdb = BioMongoDB("tdr")

    for g in genomes:
        validate_genome(g)
    print "OK"
