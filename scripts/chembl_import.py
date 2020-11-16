from tqdm import tqdm
import pandas as pd
import pymongo
import sys
import json
from collections import defaultdict

if __name__ == "__main__":
    """
    /data/databases/target/processed  
    proteinortho ../proteins.fasta ../../chembl/processed/target.fasta -cov=80
    """

    db = pymongo.MongoClient().tdr # port=27018

    props = [
        {
            "target": "protein",
            "defaultGroupOperation": "avg",
            "description": "High homology with Target in ChEMBL",
            "defaultValue": "true",
            "defaultOperation": "equal",
            "uploader": "demo",
            "type": "value",
            "options": [
                "true",
                "false"
            ],
            "name": "target_in_chembl"
        }, {
            "target": "protein",
            "defaultGroupOperation": "avg",
            "description": "High homology with protein in Reactome DB",
            "defaultValue": "true",
            "defaultOperation": "equal",
            "uploader": "demo",
            "type": "value",
            "options": [
                "true",
                "false"
            ],
            "name": "target_in_reactome"
        }, {
            "target": "protein",
            "defaultGroupOperation": "avg",
            "description": "min pChEMBL value for all avaliable activities. pChEMBL is defined as: -Log(molar IC50, XC50, EC50, AC50, Ki, Kd or Potency). For instance, a value of 5 represents an activity of 10uM and values above this are more active.",
            "defaultValue": "5",
            "defaultOperation": ">",
            "uploader": "demo",
            "type": "number",
            "options": [],
            "name": "pchembl_value"
        },
    ]
    orgs = []

    chembldata = {}
    with open("/data/databases/target/processed/chembl_data2.json") as h:
        for x in tqdm(h, total=6739):
            data = json.loads(x)
            chembldata[data["target"]["accession"]] = data

    db.proteins.update({}, {"$unset": {"search.pchembl_value": ""}, "$set": {"search.target_in_chembl": False}},
                       multi=True)
    with open("/data/databases/target/processed/chembl_blast.txt") as h:
        lines = [x.split() for x in h]
        for x in tqdm(lines):
            org = x[0].strip().split("|||")[0]
            tpp = x[0].strip().split("|||" + org)[1]
            unip = x[1].strip()
            if unip not in chembldata:
                continue
            chembl = chembldata[unip]

            update_data = {
                "chembl": {
                    "target": chembl["target"]["chemblid"],
                    "assays": []
                },
                "search.target_in_chembl": True
            }
            pval = None
            for assay in chembl["assays"]:
                dbassay = {"assay": assay["chemblid"], "activities": [], "type": assay["type"],
                           "description": assay["description"]}
                update_data["chembl"]["assays"].append(dbassay)
                for act in assay["activities"]:
                    dbact = {"compound": act["compound"]["chemblid"],
                             "comp_code":act["compound"]["key"],"ro5":act["compound"]["properties"]["num_ro5_violations"]}
                    dbassay["activities"].append(dbact)
                    if act["pchembl_value"]:
                        if pval:
                            dbact["pchembl_value"] = act["pchembl_value"]
                            if pval > act["pchembl_value"]:
                                pval = act["pchembl_value"]
                        else:
                            pval = act["pchembl_value"]
            if pval:
                update_data["search.pchembl_value"] = pval
            db.proteins.update({"organism": org, "gene": tpp}, {"$set": update_data})

            orgs.append(org)

    # db.proteins.update({}, {"$set": {"search.target_in_reactome": False, "reactome": []}}, multi=True)
    # reactome = pd.read_table("/data/databases/reactome/UniProt2ReactomeReactions.txt", sep="\t",
    #                          names=["uniprot", "reactome_id", "url", "event", "evidence", "species"])
    # react_data = {}
    # for _, r in tqdm(reactome.iterrows()):
    #     react_data[r.uniprot] = r
    #
    # with open("/data/databases/target/processed/reactome_blast_90ident.txt") as h:
    #     lines = [x.split() for x in h]
    #     for x in tqdm(lines):
    #         org = x[0].strip().split("|||")[0]
    #         orgs.append(org)
    #         tpp = x[0].strip().split("|||" )[-1]
    #         unip = x[1].strip().split("|")[1]
    #         if unip not in react_data:
    #             print (unip)
    #             continue
    #         react = react_data[unip]
    #
    #         db.proteins.update({"organism": org, "gene": tpp},
    #                            {"$set": {"search.target_in_reactome": True}})
    #         db.proteins.update({"organism": org, "gene": tpp},
    #                            {"$push":{"reactome": react.to_dict()}})
    #
    # for organism in set(orgs):
    #     for p in props:
    #         db.sequence_collection.update({"name": organism},
    #                                       {"$pull": {"druggabilityParams": {"name": p["name"]}}})
    #         db.sequence_collection.update({"name": organism},
    #                                       {"$addToSet": {"druggabilityParams": p}})
