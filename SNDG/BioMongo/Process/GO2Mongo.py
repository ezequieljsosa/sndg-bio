"""
Created on Feb 21, 2014

@author: Ezequiel
"""
import logging
import os
import copy
from _collections import defaultdict
from pymongo import MongoClient
from pymongo.database import Database
import networkx as nx
from networkx.algorithms.dag import ancestors
from goatools.obo_parser import GODag
import mongoengine
from SNDG.BioMongo.Model.Ontology import Ontology
from SNDG.BioMongo.Model.SeqColOntologyIndex import SeqColOntologyIndex
from SNDG.BioMongo.Process.KeywordIndexer import KeywordIndexer
from SNDG.Sequence.so import SO_ROOT_TERMS
from tqdm import tqdm


GO_ROOT_TERMS = ["GO:0008150", "GO:0005575", "GO:0003674"]

_log = logging.getLogger(__name__)


class GO2Mongo(object):
    """

    """

    def __init__(self, obo_file="/data/databases/go/go.obo", db="xomeq", client=None, go="ontologies",
                 go_index="col_ont_idx", ontology_name="go", slim_file=None):
        """"""
        self.obo_file = obo_file
        self.slim_file = slim_file
        self.graph = nx.DiGraph()
        self.graph_file = '/data/databases/' + ontology_name + '/' + ontology_name + '.gpickle'
        self.go_dag = None

        if isinstance(db, str):
            if not client:
                client = MongoClient()
            self.db = client[db]
        else:
            assert isinstance(db, Database)
            self.db = db

        self.col_go = self.db[go]
        self.col_go_index = self.db[go_index]
        self.ontology_name = ontology_name
        if self.ontology_name == "go":
            self.root_terms = GO_ROOT_TERMS
        else:
            self.root_terms = SO_ROOT_TERMS
        self.ki = KeywordIndexer()

    def init(self):
        _log.debug("Cargando archivo de ontologias:" + self.obo_file)
        self.go_dag = GODag(self.obo_file)
        _log.debug("Se cargo el archivo:" + self.obo_file)

        if os.path.exists(self.graph_file):
            self.graph = nx.read_gpickle(self.graph_file)
        else:
            self._build_graph()
            nx.write_gpickle(self.graph, self.graph_file)

        _log.debug("Se genero el grafo de terminos")

    def add_unknow(self):
        """
        pepe = {
            "_id" : ObjectId("591f14deaab82b7f88ef8c04"),
            "term" : "go:9999999",
            "name" : "Uknown",
            "ontology" : "go",
            "databases" : [ ],
            "description" : "Uknown",
            "keywords" : [ ],
            "parents" : [
                "go:0005575",
                "go:0003674",
                "go:0008150"
            ],
            "children" : [ ],
            "successors" : [ ],
            "subclases" : [ ],
            "successors_relationships" : [ ],
            "database" : "biological_process"
        }
        
        
        db.getCollection('ontologies').save(pepe)
        
        db.proteins.update({organism:"LactoUV",ontologies:{$ne:"go:0008150"}},{$addToSet:{ontologies:"go:9999999"}},{multi:true})
db.ontologies.update({ontology:"go",term:"go:0005575"},{$addToSet:{"subclases":"go:9999999"}})
db.ontologies.update({ontology:"go",term:"go:0005575"},{$addToSet:{"successors":"go:9999999"}})
db.ontologies.update({ontology:"go",term:"go:0005575"},{$addToSet:{"children":"go:9999999"}})


db.ontologies.update({ontology:"go",term:"go:0005575"},{$pull:{"subclases":"go:9999999"}})
db.ontologies.update({ontology:"go",term:"go:0005575"},{$pull:{"successors":"go:9999999"}})
db.ontologies.update({ontology:"go",term:"go:0005575"},{$pull:{"children":"go:9999999"}})

        
        {
    "_id" : ObjectId(""),
    "_cls" : "SeqColOntologyIndex",
    "term" : "go:9999999",
    "name" : "unknown",
    "count" : 800,
    "order" : 27,
    "keywords" : [   ],
    "ontology" : "go",
    "database" : "biological_process",
    "seq_collection_name" : "LactoUV",
    "seq_collection_id" : ObjectId("591caafebe737e774090b78d")
}
        
        
    """

    def load(self):
        self.init()
        self._load_mongo()
        _log.info("Obo %s loaded in %s collection %s and index in %s" %
                  (self.obo_file, self.db.name, self.col_go.name, self.col_go_index.name))

        _log.debug("Loading generic slim")
        self.load_slim()
        _log.info("Generic slim terms loaded")

    def load_slim(self, slim_file="/data/databases/go/goslim_generic.obo", database="generic"):
        parser = GODag(slim_file)
        for ont in parser:
            try:
                go = Ontology.objects(ontology="go", term=ont.lower()).get()
                go.databases.append(database)
                go.save()
            except Exception as ex:
                _log.error(ex)
        go = Ontology.objects(ontology="go", term="root").get()
        go.databases.append(database)
        go.save()

    def _load_mongo(self):
        root = Ontology(ontology=self.ontology_name, term="root", successors=self.root_terms,
                        children=self.root_terms)
        root.save()
        for (node, data) in self.graph.nodes_iter(data=True):  # self.graph.add_node(node, **data)
            if node == "root":
                raise Exception("...")
            else:
                successors = self.graph.successors(node)
                _ancestors = self.complete_subgraph([node])

                database = "biological_process"
                if "go:0005575" in _ancestors:
                    database = "cellular_component"
                if "go:0003674" in _ancestors:
                    database = "molecular_function"

                ont_doc = Ontology(ontology=self.ontology_name, term=node,
                                   name=data["name"],
                                   database=database,
                                   successors=self.all_successors(node, []),
                                   children=successors,
                                   description=self.go_dag.query_term(node.upper()).desc,
                                   # successors_relationships=self.successors_relationships(node),
                                   subclases=list(set(
                                       [x.lower() for x in self.go_dag.query_term(node.upper()).get_all_children()]))
                                   )
                ont_doc.keywords = self.ki.extract_keywords([ont_doc.description, ont_doc.name, ont_doc.term])
                ont_doc.save()

    def _build_graph(self):
        assert self.go_dag, "GO terms where not loaded"

        self.graph.add_node("root", name="root")

        processed = []
        for root_term in self.root_terms:  # Iterates over each root    
            root = self.go_dag.query_term(root_term)

            self.graph.add_node(root_term.lower(), name=root.name)
            self.graph.add_edge("root", root_term.lower())
            self._load_branch(root, processed)

    def _load_branch(self, term, processed):

        term_id = term.id.lower()
        if term_id in processed: return
        processed.append(term_id)
        if term.children:  # or term.relationships:
            for child in term.children:  # + [x[1] for x in term.relationships]):
                child_id = child.id.lower()

                self.graph.add_node(child_id, name=child.name)
                self.graph.add_edge(term_id, child_id)
                self._load_branch(child, processed)

    def all_successors(self, node, processed):
        successors = self.graph.successors(node)
        if successors:
            for x in successors:
                if x not in processed:
                    processed.append(x)
                    successors = list(set(successors + self.all_successors(x, processed)))
        else:
            successors = list()
        return successors

    def successors_relationships(self, node):
        term = self.go_dag.query_term(node.upper())
        return [[x.id.lower(), "is_a"] for x in term.children] + [[x[1].id.lower(), x[0]] for x in term.relationships]

    def cleanup_cellular_component_annotations(self, genome):
        for ont_doc in Ontology.objects(ontology=self.ontology_name, database="cellular_component",
                                        databases__ne="generic"):
            # self.db["proteins"].update({"organism":genome, }, {"$pull":{"ontologies":ont_doc.term, "keywords":ont_doc.term}}, multi=True)
            self.db["col_ont_idx"].remove({"ontology": "go", "seq_collection_name": genome.name, "term": ont_doc.term},
                                          multi=True)

    def complete_subgraph(self, ontologies):
        allontologies = copy.copy(ontologies)

        for ontology in ontologies:
            allontologies += self._complete_parents(ontology, [])

        return [x for x in set(allontologies) if x != "root"]

    def _complete_parents(self, ontology, walked):
        allontologies = [ontology]
        walked.append(ontology)
        if ontology in self.graph:
            for ancestor in ancestors(self.graph, ontology):
                if ancestor not in walked:
                    allontologies += self._complete_parents(ancestor, walked)
            return allontologies
        else:
            return allontologies

    #
    def pre_build_index(self, genome, annotated_collection="proteins", annotated_collection_field="ontologies",
                        drop=True):
        if drop:
            print (self.col_go_index.remove({"seq_collection_id": genome.id, "ontology": self.ontology_name}))

        ont_succ_cache = {}
        for ont_doc in tqdm(Ontology.objects(ontology=self.ontology_name).no_cache(),
                            total=Ontology.objects(ontology=self.ontology_name).count()):
            ont_succ_cache[ont_doc.term] = ont_doc.successors
            database = ""

            if hasattr(ont_doc, "database") and ont_doc.database:
                database = ont_doc.database
            #             if hasattr(ont_doc, "databases") and ont_doc.databases:
            #                 database = ont_doc.databases[0]
            order = len(ont_doc["children"])

            seq_ont_ont_idx = SeqColOntologyIndex(term=ont_doc.term.lower(), name=ont_doc.name, count=0,
                                                  seq_collection_name=genome.name, database=database,
                                                  ontology=self.ontology_name, order=order,
                                                  seq_collection_id=genome.id,
                                                  keywords=ont_doc.keywords)
            seq_ont_ont_idx.save()

        ont_count = defaultdict(lambda: 0)
        query = {"seq_collection_id": genome.id, "ontologies.0": {"$exists": True}}
        for p in tqdm(self.db[annotated_collection].find(query, {"ontologies": 1}),
                      total=self.db[annotated_collection].count(query)):
            terms = [x for x in p["ontologies"] if x.startswith("go:")]
            terms = self.complete_subgraph(terms)
            for x in terms:
                ont_count[x] += 1
            self.db[annotated_collection].update({"_id": p["_id"]
                                                  }, {"$addToSet": {annotated_collection_field: {"$each": terms}}})

        for term, count in tqdm(ont_count.items()):
            for seq_ont_ont_idx in SeqColOntologyIndex.objects(seq_collection_id=genome.id, ontology=self.ontology_name,
                                                               term=term):
                seq_ont_ont_idx.count = count
                seq_ont_ont_idx.save()

        SeqColOntologyIndex.objects(seq_collection_id=genome.id, count=0).delete()

        self.cleanup_cellular_component_annotations(genome)


if __name__ == '__main__':
    from SNDG import init_log

    init_log()
    mongoengine.connect("saureus")

    go2mongo = GO2Mongo("/data/databases/go/go-basic.obo", "saureus")
    go2mongo.db.ontologies.remove({"ontology": "go"})

    go2mongo.load()
    go2mongo.load_slim("/data/databases/go/goslim_plant.obo", "plant")
    #     go2mongo.load()
    #     go2mongo.init()
    #     genome = SeqCollection.objects(name="SaureusN315").get()
    #     go2mongo.pre_build_index(genome)
    #     go2mongo.cleanup_cellular_component_annotations(genome)
    print ("ok!")
