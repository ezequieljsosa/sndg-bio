from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
from SNDG.BioMongo.Process.Importer import index_seq_collection, build_statistics
import pymongo
from tqdm import tqdm

mdb = BioMongoDB("saureus", 27019)

## Script para aplicar el curado manual de fede
data = open("/data/organismos/ILEX_PARA2/curacion/24082018_auto.txt").read().split("#")
import re

# ecex = re.compile("^ec")
# for l in tqdm(data):
#     genes, desc, ec = [x.strip() for x in l.strip().split("\n") if x]
#     genes = genes.split("==")
#     #    try:
#     gs = [mdb.db.proteins.find_one({"organism": "ILEX_PARA", "alias": x.strip()}, {"gene": 1})["gene"][0] for x in genes if
#           x.startswith("Ilex")]
#     #    except:
#     #        print(l.strip().split("\n"))
#     ts = [x.strip() for x in genes if x.startswith("ILEX")]
#
#     for g in gs:
#         sets = {"description": desc}
#         if "Caffeine synthase" in desc:
#             num = ""
#             if len(desc.split(" ")) == 3:
#                 num = desc.split(" ")[2]
#             sets["gene"] = [g, "CS" + num]
#             sets["name"] = "CS" + num
#         print mdb.db.proteins.update({"organism": "ILEX_PARA", "gene": g},
#                            {"$pull": {"ontologies": ecex}})
#         print mdb.db.proteins.update({"organism": "ILEX_PARA", "gene": g},
#                            {"$addToSet": {"ontologies": ("ec:" + ec)},
#                             "$set": sets})
#     for t in ts:
#         print mdb.db.proteins.update({"organism": "ILEX_PARA_TRANSCRIPT", "gene": t},
#                            {"$pull": {"ontologies": ecex}})
#         print mdb.db.proteins.update({"organism": "ILEX_PARA_TRANSCRIPT", "gene": t},
#                            {"$addToSet": {"ontologies": ("ec:" + ec)},
#                             "$set": {"description": desc}})

# from SNDG.BioMongo.Process.PathwaysAnnotator import PathwaysAnnotator
# import pymongo
# pa = PathwaysAnnotator(mdb.db, "ILEX_PARA2", "/data/organismos/ILEX_PARA2/sbml/")
# pa.sbml("/data/organismos/ILEX_PARA2/sbml/small.sbml")
# pa.annotate()
#
# filter_tax = {2,
#               22,
#               29,
#               81,
#               192,
#               193,
#               194,
#               286,
#               356,
#               396,
#               429,
#               434,
#               482,
#               506,
#               543,
#               629,
#               641,
#               666,
#               724,
#               810,
#               816,
#               817,
#               838,
#               914,
#               919,
#               976,
#               1046,
#               1090,
#               1117,
#               1118,
#               1224,
#               1227,
#               1236,
#               1239,
#               1279,
#               1297,
#               1760,
#               1762,
#               1763,
#               1883,
#               2037,
#               2062,
#               2063,
#               2093,
#               28211,
#               28216,
#               28221,
#               32011,
#               32066,
#               33877,
#               33882,
#               35798,
#               39773,
#               40117,
#               40222,
#               44249,
#               55158,
#               67814,
#               68297,
#               68336,
#               72276,
#               74152,
#               82115,
#               85006,
#               85007,
#               85008,
#               85011,
#               85023,
#               85025,
#               91061,
#               91347,
#               93681,
#               135613,
#               135617,
#               135623,
#               142182,
#               171552,
#               186801,
#               186826,
#               188708,
#               191411,
#               191412,
#               200783,
#               200795,
#               200918,
#               200930,
#               200938,
#               200940,
#               201174,
#               203682,
#               203691,
#               204441,
#               204455,
#               204457,
#               267890,
#               508458,
#               544448,
#               578822,
#               1293497,
#               2157}
#
# from SNDG.BioMongo.Process.BioCyc2Mongo import BioCyc
# biocyc = BioCyc("saureus")
# biocyc.complete_pathways("ILEX_PARA2", "/data/databases/biocyc/metacyc/pathways.dat",
#                           "/data/databases/biocyc/metacyc/reactions.dat", filter_tax)
index_seq_collection(mdb.db, "ILEX_PARA", ec=True, go=False, keywords=True, organism_idx=True, pathways=True,
                     structure=False)
build_statistics(mdb.db, "ILEX_PARA")
