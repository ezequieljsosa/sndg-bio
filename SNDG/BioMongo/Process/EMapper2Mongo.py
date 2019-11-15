"""
Created on Jun 18, 2014

@author: eze
"""

import logging
import re
from tqdm import tqdm
from pymongo.database import Database
from pymongo.mongo_client import MongoClient
from bson.objectid import ObjectId

from SNDG.BioMongo.Model.Ontology import Ontology
from SNDG.BioMongo.Model.SeqColOntologyIndex import SeqColOntologyIndex
from SNDG.BioMongo.Process.KeywordIndexer import KeywordIndexer
from SNDG.BioMongo.Model.Protein import Protein

_log = logging.getLogger(__name__)

'''

'''


class EMapper2Mongosource (object):
    """
    """

    def __init__(self, db, client=None, collection="ontologies", col_index="col_ont_idx"):
        """
        Constructor
        """
        self.collection = collection

        if isinstance(db, basestring):
            if not client:
                client = MongoClient()
            self.db = client[db]
        else:
            assert isinstance(db, Database)
            self.db = db

        self.col_index = self.db[col_index]
        self.col_ec = self.db[collection]
        self.ontology_name = "ec"
        self.ki = KeywordIndexer()

    def load_enzclass(self, enzclass_file_path):

        root = Ontology(ontology=self.ontology_name, term="root", name="ec",
                        children=["ec:1.-.-.-", "ec:2.-.-.-", "ec:3.-.-.-", "ec:4.-.-.-", "ec:5.-.-.-", "ec:6.-.-.-"])
        root.save()

        with open(enzclass_file_path) as enzclass_handle:
            for line in enzclass_handle:
                if re.match(r'^[1-6][.]', line):
                    name = line.split(".-")[-1].strip()
                    term = "ec:" + line.replace(name, "").replace(" ", "").strip()

                    ont_doc = Ontology(ontology=self.ontology_name, term=term, name=name)
                    ont_doc.keywords = self.ki.extract_keywords(ont_doc.name) + [ont_doc.term]
                    ont_doc.save()

    def load_enzdata(self, enzdata_file_path):
        ont_doc = None
        with open(enzdata_file_path) as enzclass_handle:
            for line in enzclass_handle:
                if line.startswith("DE"):
                    ont_doc.name = line.split("DE")[1].strip()
                    ont_doc.keywords = self.ki.extract_keywords([ont_doc.description, ont_doc.name]) + [ont_doc.term]
                    ont_doc.save()
                elif line.startswith("ID"):
                    term = "ec:" + line.split("ID")[1].strip()
                    ont_doc = Ontology(ontology=self.ontology_name, term=term)

    def load_children(self):
        for ont_doc in Ontology.objects(ontology="ec"):
            if ont_doc.term == "root":
                continue
            ont_doc.children = [x["term"] for x in Ontology.objects(
                ontology="ec", term__istartswith=ont_doc.term.replace(".-", "")
            ) if (x.term != ont_doc.term) and (x.term.count("-") == (ont_doc.term.count("-") - 1))]
            ont_doc.__repr__()
            ont_doc.save()

    def load_priam_hits(self, seq_collection_name, path_genomeEnzymes):

        for line in open(path_genomeEnzymes):
            arr_line = line.split("\t")
            ec = "ec:" + arr_line[0]
            prot_id = arr_line[1].split(" ")[0]

            for prot in Protein.objects(organism=seq_collection_name, alias=prot_id):
                prot.ontologies.append(ec)
                prot.keywords.append(ec)
                prot.save()

    def pre_build_index(self, genome, annotated_collection="proteins", annotated_collection_field="ontologies",
                        drop=False):
        self.col_index.remove({"ontology": "ec", "seq_collection_id": genome.id})

        _log.debug("creating empty ecs idxs")
        ec_count = Ontology.objects(ontology="ec").count()
        for ont_doc in tqdm(Ontology.objects(ontology="ec").no_cache(), total=ec_count):
            order = 9999999
            try:
                order = order / int(
                    ont_doc.term.lower().replace("ec:", "").replace(".", "").replace("-", "").replace("n", ""))
            except:
                pass

            seq_col_ont_idx = SeqColOntologyIndex(seq_collection_id=genome.id, term=ont_doc.term.lower(),
                                                  seq_collection_name=genome.name, name=ont_doc.name,
                                                  ontology="ec", order=order,
                                                  count=0, keywords=ont_doc.keywords)
            seq_col_ont_idx.save()



        _log.debug("initializign idx ecs")
        terms_count = SeqColOntologyIndex.objects(ontology="ec", seq_collection_id=genome.id).count()
        for ont_doc in tqdm(SeqColOntologyIndex.objects(ontology="ec", seq_collection_id=genome.id).order_by(
                "-term").no_cache().timeout(False), total=terms_count):
            if "-" in ont_doc.term:
                str_term = ont_doc.term.replace(".", "\.").replace("-", ".+")
            else:
                str_term = ont_doc.term

            count = self.db[annotated_collection].count({"seq_collection_id": genome.id,
                                                         annotated_collection_field: {"$regex": '^' + str_term,
                                                                                      "$options": "-i"}})
            if count:
                ont_doc.count = count
                ont_doc.save()
            else:
                ont_doc.delete()

        regx = re.compile("^ec:", re.IGNORECASE)
        self.db.col_ont_idx.insert(
            {
                "_id": ObjectId(),
                "_cls": "SeqColOntologyIndex",
                "term": "root",
                "name": "root",
                "count": self.db.proteins.count({"organism": genome.name, "ontologies": regx}),
                "order": 9999999,
                "keywords": [
                ],
                "ontology": "ec",
                "seq_collection_name": genome.name,
                "seq_collection_id": genome.id
            })
        _log.debug("ecs idxs done")

    def _process_children(self, seq_collection_id, counts, str_term):
        for successor in self._successors(str_term, counts):
            suc_count = self.col_index.find_one({"seq_collection_id": seq_collection_id, "term": successor})[
                "count"]  # counts[successor.lower()] if successor.lower() in counts else 0
            if suc_count:
                self.col_index.update({"seq_collection_id": seq_collection_id, "term": str_term},
                                      {"$inc": {"count": suc_count}})

    def _successors(self, term, counts):
        try:
            ec_term = Ontology.objects(ontology="ec", term=term).get()
            return ec_term.children  # [ x for x in counts.keys() if re.match(term.split("-")[0], x)   ]
        except:
            return []
