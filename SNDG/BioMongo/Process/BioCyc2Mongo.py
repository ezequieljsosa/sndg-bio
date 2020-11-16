'''
Created on Jun 18, 2014

@author: eze
'''

import logging
import re

import mongoengine
from SNDG.BioMongo.Model.Pathway import ChokepointType
from mongoengine.fields import ListField, StringField
from pymongo.database import Database
from pymongo.mongo_client import MongoClient

from SNDG.BioMongo.Model.Ontology import Ontology
from SNDG.BioMongo.Model.Pathway import PathwaySumary
from SNDG.BioMongo.Model.Protein import Protein
from SNDG.BioMongo.Model.SeqColDruggabilityParam import SeqColDruggabilityParam
from SNDG.BioMongo.Model.SeqColOntologyIndex import SeqColOntologyIndex
from SNDG.BioMongo.Model.Sequence import ProteinDruggabilitySearch
from SNDG.BioMongo.Process.KeywordIndexer import KeywordIndexer

from tqdm import tqdm

_log = logging.getLogger(__name__)

from SNDG.BioMongo.Model.SeqColDruggabilityParam import SeqColDruggabilityParamTypes

"""
Para ver las jerarquias entre las clases hay que usar el archivo classes.dat
"""


class BioCyc(object):
    '''
     
    '''

    protein_pathway_search_params = [
        ("centrality",
         """Shortest-path betweenness centrality (normalized) for reactions. 
         In the used graph the nodes are the reactions and the edges the metabolites conecting them.
         When centrality >= 0.1 the reaction is considered highly central
         """,
         "pathways",
         SeqColDruggabilityParamTypes.number,
         None,
         "max", ">", "0.1"),
        ("chokepoint", "If the protein catalyzes a reaction that is a metabolic chokepoint", "pathways",
         SeqColDruggabilityParamTypes.value, ["true", "false"],
         "avg", "equal", "true"),
        ("chokepoint_type", "metabolic chokepoint type, can be of consume, production or double", "pathways",
         SeqColDruggabilityParamTypes.value,
         [str(x.value) for x in ChokepointType],
         "avg", "equal", "true")]
    '''
    name,description,target,_type,options,defaultGroupOperation,defaultOperation,defaultValue
    '''

    pathways_search_params = [
        ("reactions", "Number of reactions in the pathway", "pathway", SeqColDruggabilityParamTypes.number, None, "avg",
         ">", 0),
        ("norm_reactions", "Number of reactions normalized by the highest pathway (with more reactions)", "pathway",
         SeqColDruggabilityParamTypes.number, None, "avg", ">", 0),
        ("reactions_with_gene", "Number of reactions in the pathway with at least one known protein that catalize them",
         "pathway", SeqColDruggabilityParamTypes.number, None, "avg", ">", 0),
        ("completeness", "Proportion of reactions in the pathway with at least one known protein that catalize them",
         "pathway",
         SeqColDruggabilityParamTypes.number, None, "max", ">", 0.8),
        (
            "chokepoints", "Number of chokepoints reactions in the pathway", "pathway",
            SeqColDruggabilityParamTypes.number,
            None, "avg", ">", 0),
        ("norm_chokepoint", "Chokepoint reactions/reactions in pathway ratio", "pathway",
         SeqColDruggabilityParamTypes.number, None, "avg", ">", 0),
        ("druggable", "The pathway has at least one druggable protein", "pathway", SeqColDruggabilityParamTypes.value,
         ["Yes", "No"],
         "avg", "equal", "Yes"),
        ("max_centrality",
         "Maximum betweenness centrality of all the reactions in the pathway, normalized by the reaction with max betweenes centrality in the whole graph",
         "pathway", SeqColDruggabilityParamTypes.number, None,
         "avg", ">", "0.2"),
        #         ("human_pw_offtarget", "Average reactions count without any protein with an humman offtarget < 0.6, normalized by the number of reactions with genes", "pathway",
        #          SeqColDruggabilityParamTypes.number, None,      "avg", "equal", 0)
    ]

    def __init__(self, db, client=None, collection="ontologies"):
        '''
        Constructor
        '''
        self.collection = collection

        if isinstance(db, str):
            if not client:
                client = MongoClient()
            self.db = client[db]
        else:
            assert isinstance(db, Database)
            self.db = db
        self.user = "demo"
        self.col_biocyc = self.db[collection]
        self.ontology_name = "biocyc"
        self.user = "demo"
        self.ki = KeywordIndexer()
        self.pathways = {}
        self.react_ont_idx_dict = {}

    def _protein_keywords(self, protein):
        keywords = []
        for reaction in protein.reactions:
            keywords.append(reaction.name)
            if reaction.name not in self.react_ont_idx_dict:
                seq_ont_ont_idx = SeqColOntologyIndex(term=reaction.name, seq_collection_name=protein.organism,
                                                      count=1, seq_collection_id=protein.seq_collection_id,
                                                      ontology=self.ontology_name + "_reac",
                                                      keywords=[reaction.name])
                self.react_ont_idx_dict[reaction.name] = seq_ont_ont_idx
            else:
                seq_ont_ont_idx = self.react_ont_idx_dict[reaction.name]
                seq_ont_ont_idx.count = seq_ont_ont_idx.count + 1
            keywords = keywords + reaction.pathways
            for pathway in reaction.pathways:
                if pathway not in self.pathways:
                    self.pathways[pathway] = []
                self.pathways[pathway].append(protein.id)

            for specie in reaction.substrates + reaction.products:
                keywords.append(specie.name)
                seq_ont_ont_idx.keywords.append(specie.name)
                keywords = keywords + specie.producers + specie.consumers
        return keywords

    def update_ont_org_idx(self):
        for term, seq_ont_ont_idx in self.react_ont_idx_dict.items():
            reac = Ontology.objects(ontology=self.ontology_name + "_reac", term=term.lower())
            if len(reac):
                reac = reac.first()
                seq_ont_ont_idx.name = reac.name
                seq_ont_ont_idx.keywords = list(set(seq_ont_ont_idx.keywords + reac.keywords))
                seq_ont_ont_idx.save()

    def _add_drugability_props_to_protein(self, protein):
        keywords = []
        if (not hasattr(protein, "search")) or (protein.search == None):
            protein.search = ProteinDruggabilitySearch()
        protein.search.chokepoint = False
        props = [prop for prop in protein.properties if (prop._type == "pathways") and (prop.property == "chokepoint")]
        if props:
            prop = props[0]
            protein.search.chokepoint_type = prop.type
            protein.search.chokepoint = True
            keywords.append("chokepoint")

        props = [prop for prop in protein.properties if (prop._type == "pathways") and (prop.property == "centrality")]
        if props:
            prop = props[0]
            protein.search.centrality = prop.value
        return keywords

    def _add_drugability_props_to_genome(self, genome):

        for name, description, target, _type, options, dgo, do, dv in BioCyc.protein_pathway_search_params + BioCyc.pathways_search_params:
            if not genome.has_druggability_param(name):
                dp = SeqColDruggabilityParam(name=name, description=description, target=target,
                                             type=_type, uploader=self.user)
                if options:
                    dp.options = options

                dp.defaultGroupOperation = dgo
                dp.defaultOperation = do
                dp.defaultValue = dv

                genome.druggabilityParams.append(dp)

    def _process_proteins(self, genome):
        total = Protein.objects(organism=genome.name, reactions__0__exists=True).count()
        iterprot = tqdm(Protein.objects(organism=genome.name, reactions__0__exists=True).no_cache().timeout(
            False),total=total)
        for protein in iterprot:  # @UndefinedVariable
            keywords1 = self._protein_keywords(protein)
            keywords2 = self._add_drugability_props_to_protein(protein)
            protein.keywords = list(
                set([x.strip().lower() for x in keywords1 + keywords2 + protein.keywords if x.strip()]))
            protein.save()
        self.update_ont_org_idx()

    def _genome_summary(self, genome):
        pathways_count = {x: len(set(y)) for x, y in self.pathways.items()}



        SeqColOntologyIndex.objects(seq_collection_name=genome.name, ontology=self.ontology_name + "_pw").delete()

        for x in genome.pathways :
            seq_ont_ont_idx = SeqColOntologyIndex(term=x["term"], name=x["name"],
                                                  seq_collection_name=genome.name,
                                                  ontology=self.ontology_name + "_pw",
                                                  keywords=self.ki.extract_keywords(x["name"]) + [x["term"]],
                                                  count=x["count"], seq_collection_id=genome.id)
            seq_ont_ont_idx.save()

    def pre_build_index(self, genome):
        self._process_proteins(genome)
        self._genome_summary(genome)
        self._add_drugability_props_to_genome(genome)
        for x in genome.druggabilityParams:
            x.uploader = "demo"
        genome.save()

    def load_pathways(self, pathways_file, database):
        with open(pathways_file) as pathways_handle:
            for x in pathways_handle.readlines():
                if not (x.strip().startswith("#") or x.strip().startswith("UNIQUE-ID")):
                    line = re.sub(r'\t+', '\t', x)
                    term = line.split("\t")[0].strip().lower()
                    name = line.split("\t")[1].strip()
                    ont_doc = Ontology(term=term,
                                       name=name,
                                       ontology=self.ontology_name + "_pw",
                                       keywords=self.ki.extract_keywords(name) + [term])
                    ont_doc.databases.append(database)
                    ont_doc.save()

    def complete_pathways(self, genome, pathways_file, reactions_file, filter_tax=None):
        """

        :param genome:
        :param pathways_file:
        :param reactions_file:
        :param filter_tax: list of ncbi_taxon_id (only number)
        :return:
        """
        if filter_tax:
            filter_tax = set([str(x) for x in filter_tax])
        pathways = []
        for p in self.db.proteins.find({"organism": genome, "reactions.pathways.0": {"$exists": 1}}, {"reactions": 1}):
            for rs in p["reactions"]:
                pathways = pathways + rs["pathways"]

        pathways = {x: 1 for x in pathways}

        gecs = []
        regx = re.compile("^ec:")
        for p in self.db.proteins.find({"organism": genome, "ontologies": regx}, {"ontologies": 1}):
            gecs = gecs + [x for x in p["ontologies"] if x.startswith("ec:")]
        gecs = set(gecs)
        react_ec = {}

        with open(reactions_file) as h:
            lines = [x for x in h.readlines() if not x.startswith("#")]
            records = re.split("//\n", "\n".join(lines))
            for record in tqdm(records):
                if not record.strip():
                    continue
                ec = None
                uid = None
                for str_record in [y for y in record.split("\n") if y]:
                    if str_record.strip() and len(str_record.strip()) > 3:

                        if len(str_record.split(" - ")) > 1:

                            field = str_record.split(" - ")[0].strip()
                            try:
                                value = str_record.split(" - ")[1].strip().decode("utf-8")
                            except UnicodeDecodeError:
                                continue

                            if field == "UNIQUE-ID":
                                uid = value.lower()
                            elif (field == "EC-NUMBER"):
                                ec = value
                if ec and uid:
                    react_ec[uid] = ec.replace("-", ":").lower()

        with open(pathways_file) as h:
            lines = [x for x in h.readlines() if not x.startswith("#")]
            records = re.split("//\n", "\n".join(lines))
            for record in tqdm(records):
                if not record.strip():
                    continue
                reactions = []
                uid = None
                taxs = []
                name = None
                for str_record in [y for y in record.split("\n") if y]:
                    if str_record.strip() and len(str_record.strip()) > 3:
                        if len(str_record.split(" - ")) > 1:
                            field = str_record.split(" - ")[0].strip()
                            try:
                                value = str_record.split(" - ")[1].strip().decode("utf-8")
                            except UnicodeDecodeError:
                                continue

                            if field == "UNIQUE-ID":
                                uid = value.lower()

                            elif (field == "REACTION-LIST"):
                                reactions.append(value)
                            elif (field == "TAXONOMIC-RANGE"):
                                taxs.append(value.split("-")[1])
                            elif (field == "COMMON-NAME"):
                                name = value
                ecs = []
                if uid and reactions and uid.upper() not in pathways:
                    for r in reactions:
                        if r.lower() in react_ec:
                            ecs.append(react_ec[r.lower()])
                    ecs = set(ecs)
                    taxs = set(taxs)
                    if ecs and (len(reactions) > 3) and ((len(gecs & ecs) * 1.0 / len(reactions)) > 0.5):
                        if not filter_tax or len(taxs - filter_tax):

                            for r in reactions:
                                for p in self.db.proteins.find({"organism": genome, "reactions.name": r}):
                                    reac = [x for x in p["reactions"] if x["name"] == r][0]
                                    reac["pathways"].append(uid.upper())
                                    self.db.proteins.save(p)

    def load_dat(self, reactions_file, database, postfix):
        with open(reactions_file) as reactions_handle:
            lines = [x for x in reactions_handle.readlines() if not x.startswith("#")]
            records = re.split("//\n", "\n".join(lines))
            for record in records:
                if not record.strip():
                    continue

                ont_doc = Ontology(ontology=self.ontology_name + postfix)
                ont_doc.databases.append(database)
                reaction_types = []
                ec = None
                for str_record in [y for y in record.split("\n") if y]:
                    if str_record.strip() and len(str_record.strip()) > 3:

                        if len(str_record.split(" - ")) > 1:

                            field = str_record.split(" - ")[0].strip()
                            try:
                                value = str_record.split(" - ")[1].strip().decode("utf-8")
                            except UnicodeDecodeError:
                                continue

                            if field == "UNIQUE-ID":
                                ont_doc.term = value.lower()
                            elif field == "TYPES":
                                reaction_types.append(value)
                            elif field == "IN-PATHWAY":
                                ont_doc.parents.append(value)
                            elif field == "COMMON-NAME":
                                ont_doc.name = value
                            elif (field == "COMMENT") and (not ont_doc.name):
                                ont_doc.description = value
                            elif (field == "EC-NUMBER") and (not ont_doc.name):
                                ec = value

                if not ont_doc.description:
                    ont_doc.description = "|".join(reaction_types)
                if not ont_doc.name:
                    if ec:
                        ont_doc.name = ec
                    else:
                        ont_doc.name = ont_doc.term
                ont_doc.keywords = self.ki.extract_keywords(ont_doc.name) + [ont_doc.term]
                ont_doc.types = reaction_types
                if ec:
                    ont_doc.keywords.append(ec)
                if not ont_doc.term:
                    print (record)
                else:
                    ont_doc.save()

    def load_reactions(self, reactions_file, database):
        self.load_dat(reactions_file, database, "_reac")

    def load_compounds(self, reactions_file, database):
        self.load_dat(reactions_file, database, "_comp")

    def load_pathways_meta(self, pathways_file):
        with open(pathways_file) as pathways_handle:
            lines = [x for x in pathways_handle.readlines() if not x.startswith("#")]
            records = re.split("//\n", "\n".join(lines))
            for record in records:
                if not record.strip():
                    continue

                ont_doc = BiocycPWOnt(ontology=self.ontology_name + "_meta")
                reaction_types = []

                for str_record in [y for y in record.split("\n") if y]:
                    if str_record.strip() and len(str_record.strip()) > 3:

                        if len(str_record.split(" - ")) > 1:

                            field = str_record.split(" - ")[0].strip()
                            try:
                                value = str_record.split(" - ")[1].strip().decode("utf-8")
                            except UnicodeDecodeError:
                                continue

                            if field == "UNIQUE-ID":
                                ont_doc.term = value.lower()
                            elif field == "TYPES":
                                reaction_types.append(value)
                            elif field == "IN-PATHWAY":
                                ont_doc.parents.append(value)
                            elif field == "COMMON-NAME":
                                ont_doc.name = value
                            elif (field == "COMMENT"):
                                ont_doc.description = value
                            elif (field == "PATHWAY-LINKS"):
                                if "PWY-" in value:
                                    for x in re.findall("[a-zA-Z\-]+", value):
                                        if "PWY" in x:
                                            ont_doc.pwlinks.append(x)
                                else:
                                    value = re.findall("[a-zA-Z\-]+", value)[-1]
                                    ont_doc.pwlinks.append(value)

                            elif (field == "TAXONOMIC-RANGE"):
                                try:
                                    ont_doc.taxs.append(value.split("TAX-")[1])
                                except:
                                    pass
                            elif (field == "REACTION-LIST"):
                                ont_doc.reactions.append(value)

                if not ont_doc.description:
                    ont_doc.description = "|".join(reaction_types)
                if not ont_doc.name:
                    ont_doc.name = ont_doc.term
                ont_doc.keywords = self.ki.extract_keywords(ont_doc.name) + [ont_doc.term]
                ont_doc.types = reaction_types
                if not ont_doc.term:
                    print (record)
                else:
                    ont_doc.save()


class BiocycPWOnt(Ontology):
    '''
    
    '''

    reactions = ListField(StringField(), default=[])
    taxs = ListField(StringField(), default=[])
    pwlinks = ListField(StringField(), default=[])
    types = ListField(StringField(), default=[])


if __name__ == '__main__':
    bacs = ['666', '85007', '82115', '2093', '482', '914', '203682', '40222', '29', '543', '85011', '816', '28211',
            '1090', '33882', '286', '85025', '200938', '2062', '2063', '85006', '976', '203691', '194', '193', '192',
            '32011', '396', '85008', '93681', '919', '91347', '1236', '135613', '1239', '200795', '40117', '188708',
            '810', '28221', '135617', '85023', '22', '200940', '641', '434', '1279', '1227', '44249', '1224', '629',
            '72276', '1883', '2', '506', '55158', '191411', '74152', '1293497', '171552', '33877', '135623', '204457',
            '32066', '186826', '191412', '544448', '724', '200783', '578822', '68336', '91061', '68297', '39773',
            '200918', '267890', '81', '204441', '28216', '201174', '2037', '67814', '1117', '1118', '817', '1297',
            '186801', '142182', '200930', '429', '35798', '1046', '204455', '356', '508458', '1762', '1763', '1760',
            '838']
    arc = ["2157"]
    mongoengine.connect("saureus")
    bcyc = BioCyc("saureus")

    bcyc.db.ontologies.remove({"ontology": "biocyc_meta"})
    bcyc.load_pathways_meta("/data/databases/biocyc/metacyc/pathways.dat")
    #     bcyc.db.ontologies.remove({"ontology":"biocyc_comp"})
    #     bcyc.load_compounds("/data/databases/biocyc/metacyc/compounds.dat", "metacyc")
    #     bcyc.load_compounds("/data/databases/biocyc/ecocyc/compounds.dat", "ecocyc")

    print ("OK")
