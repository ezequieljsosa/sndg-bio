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
db = MongoClient().saureus

import re

regx_go = re.compile("^go:")  # , re.IGNORECASE)
regx_ec = re.compile("^ec:")

to_correct = []

def validate_prots(g):
    gos = db.proteins.count({"organism": g.name, "ontologies": regx_go})
    if gos < 500 and db.proteins.count({"organism": g.name}) > 100:
        print "%s tiene pocas proteinas con go anotados: %i" % (g.name, gos)
        to_correct.append(g.name)

    ecs = db.proteins.count({"organism": g.name, "ontologies": regx_ec})
    if ecs < 500 and db.proteins.count({"organism": g.name}) > 100:
        print "%s tiene pocas proteinas con ec anotados: %i" % (g.name, ecs)


def validate_genome(g):

    validate_prots(g)
    for x in ["ec", "go"]:
        if db.col_ont_idx.count({"ontology": x, "seq_collection_name": g.name}) == 0:
            print  g.name + " sin indice " + x



if __name__ == '__main__':
    BioMongoDB("saureus")
    genomes = list(Genome.objects(auth=BioMongoDB.demo_id))

    assert 100 < len(genomes), len(genomes)

    no_stats = db.sequence_collection.count({"statistics.0": {"$exists": False}})
    if no_stats:
        print "there are %i genomes with no stats!!" % no_stats
    for g in genomes:
        validate_genome(g)

    print "-------------"
    print to_correct
    print "OK"
