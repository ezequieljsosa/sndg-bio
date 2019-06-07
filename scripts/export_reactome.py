from datetime import datetime
from tqdm import tqdm
import pandas as pd
import pymongo
import sys
import json
from collections import defaultdict

db = pymongo.MongoClient(port=27018).tdr
dateTimeObj = datetime.now()
timestampStr = dateTimeObj.strftime("%Y_%m_%d")

with open('/data/databases/target/reactome_map_' + timestampStr + '.txt', "w") as h:
    h.write("target_url,genome,target_genes,reactome_id,uniprot\n")
    for p in tqdm(db.proteins.find({"search.target_in_reactome":True}, {"organism": 1, "gene": 1, "reactome": 1}),
                  total=db.proteins.count({"search.target_in_reactome":True})):
        for r in p["reactome"]:
            h.write(",".join(
                ["http://target.sbg.qb.fcen.uba.ar/patho/protein/" + str(p["_id"]), str(p["organism"]), "|".join(p["gene"]) ,
                 str(r["reactome_id"]), str(r["uniprot"])]
            ) + "\n")
