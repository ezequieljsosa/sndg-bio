"""
Created on May 8, 2015

@author: eze
"""
import re

from pymongo.mongo_client import MongoClient

from SNDG.BioMongo.Model.Ontology import Ontology
from SNDG.BioMongo.Process.KeywordIndexer import KeywordIndexer


class COG2Mongo(object):

    def __init__(self):

        self.ki = KeywordIndexer()

    def load_whog(self, whog_path):
        with open(whog_path, "r") as handle:
            groups = handle.read().split("_______")
        for group in groups:
            glines = [g.strip() for g in group.split("\n") if g.strip()]
            if glines:
                line = glines[0].lower()
                parent = line.lower().strip().split(" ")[0]
                term = line.lower().strip().split(" ")[1]
                name = " ".join(line.lower().strip().split(" ")[2:])

                ont_doc = Ontology(term=term, name=name, parent=parent,
                                   ontology="cog")
                keywords = self.ki.extract_keywords(line)

                if len(parent) > 3:
                    for x in parent[1:-1]:
                        parent_ont_doc = Ontology.objects(term='[' + x + ']').get()
                        keywords = list(set(parent_ont_doc.keywords + keywords))
                        parent_ont_doc.children.append(term)
                        parent_ont_doc.save()
                else:
                    parent_ont_doc = Ontology.objects(term=parent).get()
                    parent_ont_doc.children.append(term)
                    parent_ont_doc.save()
                    keywords = list(set(parent_ont_doc.keywords + keywords))

                ont_doc.keywords = keywords

                ont_doc.save()

    def load_fun(self, fun_path):
        with open(fun_path, "r") as handle:
            string = handle.read()
            groups = re.split(re.compile("\r\n\r\n"), string)
        for group in groups:
            glines = [g.strip() for g in group.split("\n") if g.strip()]

            term = glines[0].strip().lower()

            keywords = self.ki.extract_keywords(term)
            parent_ont_doc = Ontology(term=term, name=term,
                                      keywords=keywords, ontology="cog")
            parent_ont_doc.save()

            for line in glines[1:]:
                term = line.lower().strip().split(" ")[0]
                name = " ".join(line.lower().strip().split(" ")[1:])
                keywords = list(set(parent_ont_doc.keywords + self.ki.extract_keywords(line)))
                ont_doc = Ontology(term=term, name=name, parent=parent_ont_doc.term,
                                   keywords=keywords, ontology="cog")
                parent_ont_doc.children.append(term)
                ont_doc.save()

    #     def load_go(self,whoggo_path):
    #         with open(whoggo_path,"r") as handle:
    #             for line in handle.readlines():
    #                 if line.strip() and not line.startswith("!") and not line.endswith("GO:."):
    #                     go = line.strip().split(";")[-1].strip()
    #                     cog = line.strip().split(">")[0].replace("COG:","")[1:].strip()
    #                     keywords = cog

    def create_ontology(self, ontology_db):
        ontology_db.remove({"ontology": "cog"})
        for ontology in self.cog:
            keywords = self.ki.extract_keywords(ontology)
            ont_doc = Ontology(term=ontology,
                               keywords=keywords, ontology="cog")
            ont_doc.save()


if __name__ == '__main__':
    cog = COG2Mongo()
    cog.load_whog("/data/cog/whog")
    cog.create_ontology(MongoClient().test_database.ontologies)
