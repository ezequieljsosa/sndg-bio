from SNDG.BioMongo.Process.Importer import load_pathways,build_statistics
from SNDG.BioMongo.Process.BioCyc2Mongo import BioCyc
from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
from SNDG.BioMongo.Model.SeqCollection import SeqCollection
import pymongo

port = 27018
mdb = BioMongoDB("tdr",port=port)
db = pymongo.MongoClient(port=port).pdb

load_pathways("cruzi", "/data/organismos/cruzi/pathways/pathways-sm.sbml", mdb.db,
             "/data/organismos/cruzi/pathways/", filter_file="allfilters_con_c.dat")
biocyc = BioCyc(mdb.db)
biocyc.user = BioMongoDB.demo
biocyc.pre_build_index(SeqCollection.objects(name="cruzi").get())
build_statistics(mdb.db,"cruzi")