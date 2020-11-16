'''
Created on Sep 23, 2016

@author: eze
'''
import itertools
import logging
import os
import pickle
from _collections import defaultdict
import networkx as nx
from tqdm import tqdm
import json
from functools import reduce
from mongoengine.errors import DoesNotExist, MultipleObjectsReturned
from networkx.algorithms.centrality.betweenness import betweenness_centrality
from networkx.algorithms.components.connected import connected_components

from SNDG.BioMongo.Model.Pathway import Reaction, ChemSpecie, PathwaySumary, ChokepointType
from SNDG.BioMongo.Model.Protein import Protein
from SNDG.BioMongo.Model.Sequence import BioProperty
from SNDG.Network.SBMLProcessor import SBMLProcessor

_log = logging.getLogger(__name__)


class PathwaysAnnotator(object):
    '''

    '''

    def __init__(self, db, organism, work_dir):
        '''
        '''
        self.sbmlprocessor = SBMLProcessor()
        self.db = db
        self.prefered_biocyc = None
        self.organism = organism
        self.work_dir = work_dir
        self.unmapped_genes = []
        self.pathways = []

    def sbml(self, sbml_file):
        self.sbmlprocessor.sbml_file_path = sbml_file
        return self

    def species_filter(self, species_filter_file):
        self.sbmlprocessor.filter_filename = species_filter_file
        return self

    def chokepoints(self):
        metabolites_in = defaultdict(list)
        metabolites_out = defaultdict(list)
        for n, data in self.sbmlprocessor.graph.nodes(data=True):
            for m in data["products"]:
                if m not in self.sbmlprocessor.filter_list:
                    metabolites_out[m].append(n)
            for m in data["reactants"]:
                if m not in self.sbmlprocessor.filter_list:
                    metabolites_in[m].append(n)

        def extract_genes(met_dict):
            return set(
                reduce(list.__add__, [list(itertools.product(self.sbmlprocessor.genes[reactions[0]], [metabolite]))
                                      for metabolite, reactions in met_dict.items()
                                      if len(reactions) == 1]))

        def extract_reactions(met_dict):
            return set([reactions[0]
                        for _, reactions in met_dict.items()
                        if len(reactions) == 1
                        and len(self.sbmlprocessor.genes[reactions[0]]) == 1])

        self._chokepoints = list(
            set(list(extract_reactions(metabolites_in)) + list(extract_reactions(metabolites_out))))
        self.consuming_chokepoints, self.production_chokepoints = extract_genes(metabolites_in), extract_genes(
            metabolites_out)
        self.double_chokepoints = self.consuming_chokepoints & self.production_chokepoints

    def extract_genes_from_notes(self, fn_extract):
        self.sbmlprocessor.fn_extract_genes = fn_extract
        return self

    #     def add_pathways_info_to_gene(self, reaction, protein):
    #         if reaction in self.sbmlprocessor.pathways:
    #             if hasattr(protein, "reactions") :
    #                 protein.reactions = []
    #             for pathway_raw in self.sbmlprocessor.pathways[reaction]:
    #                 pathway = pathway_raw.replace(".","_")
    #                 gene_dict["pathways"][pathway] = {}

    def add_reaction_info_to_gene(self, reaction_raw, protein):
        '''
        "protein" is the  enzime that catalizes the "reaction_raw" 
        '''
        reaction_name = reaction_raw.replace(".", "_")
        reaction = protein.reaction(name=reaction_name)
        if not reaction:
            raw_pathways = self.sbmlprocessor.pathways[reaction_raw]
            pathways = [pathway_raw.replace(".", "_") for pathway_raw in raw_pathways]
            reaction = Reaction(name=reaction_name, substrates=[], products=[], pathways=pathways)
            protein.reactions.append(reaction)

        connected_reactions = list(nx.all_neighbors(self.sbmlprocessor.graph, reaction_raw))
        if connected_reactions:
            connected_reactions_with_genes = set([neighbor_reaction for neighbor_reaction in connected_reactions if
                                                  neighbor_reaction in self.sbmlprocessor.genes])
            if connected_reactions_with_genes:
                for neighbor_reaction_raw in connected_reactions_with_genes:
                    neighbor_reaction = neighbor_reaction_raw.replace(".", "_")
                    reaction_genes = self.sbmlprocessor.genes[neighbor_reaction_raw]
                    if reaction_genes:
                        if reaction_raw in self.sbmlprocessor.pathways:
                            # Add product and genes to reaction
                            if reaction_raw in self.sbmlprocessor.graph \
                                    and neighbor_reaction_raw in self.sbmlprocessor.graph[reaction_raw]:
                                # and "linkedby" in self.sbmlprocessor.graph[reaction_raw][neighbor_reaction_raw]:
                                specie_name = self.sbmlprocessor.graph[reaction_raw][neighbor_reaction_raw]["linkedby"]
                                specie = reaction.product(name=specie_name)

                                if specie:
                                    specie.producers = list(set(specie.producers + reaction_genes))
                                else:
                                    specie = ChemSpecie(name=specie_name, producers=list(reaction_genes))
                                    reaction.products.append(specie)

                                    # Add sbutrates and genes to reaction
                            if neighbor_reaction_raw in self.sbmlprocessor.graph \
                                    and reaction_raw in self.sbmlprocessor.graph[neighbor_reaction_raw]:
                                # and "linkedby" in self.sbmlprocessor.graph[neighbor_reaction_raw][reaction_raw]:
                                specie_name = self.sbmlprocessor.graph[neighbor_reaction_raw][reaction_raw]["linkedby"]
                                specie = reaction.substrate(name=specie_name)

                                if specie:
                                    specie.consumers = list(set(specie.consumers + reaction_genes))
                                else:
                                    specie = ChemSpecie(name=specie_name, consumers=list(reaction_genes))
                                    reaction.substrates.append(specie)

    def add_properties_to_gene(self, protein, gene_name):
        gene = gene_name
        choke_type = None
        double_cp_genes = [x[0] for x in self.double_chokepoints]
        consuming_cp_genes = [x[0] for x in self.consuming_chokepoints]
        productorion_cp_genes = [x[0] for x in self.production_chokepoints]

        def add_choke(choke_type, met):
            prop = BioProperty(_type="pathways", property="chokepoint", type=choke_type, metabolites=met)
            if not protein.has_property(prop):
                protein.properties.append(prop)

        if gene in double_cp_genes:
            met = [m for cpgene, m in self.double_chokepoints if cpgene == gene]
            choke_type = ChokepointType.double.value  # @UndefinedVariable        
            add_choke(choke_type, met)
        elif gene in consuming_cp_genes:
            met = [m for cpgene, m in self.consuming_chokepoints if cpgene == gene]
            choke_type = ChokepointType.consuming.value  # @UndefinedVariable
            add_choke(choke_type, met)
        elif gene in productorion_cp_genes:
            met = [m for cpgene, m in self.production_chokepoints if cpgene == gene]
            choke_type = ChokepointType.production.value  # @UndefinedVariable
            add_choke(choke_type, met)

        max_centrality = max(self.centrality.values())
        centrality_data = [self.centrality[str(r).replace("_", ".")] for r in self.gene_reaction[gene_name] if
                           str(r).replace("_", ".") in self.centrality]
        if centrality_data and max_centrality:
            centrality = max(centrality_data) / max_centrality
            prop = BioProperty(_type="pathways", property="centrality", value=centrality)
            if not protein.has_property(prop):
                protein.properties.append(prop)

    def build_network(self):

        if not os.path.exists(self.work_dir + "/" + self.sbmlprocessor.filter_filename):
            self.sbmlprocessor.init()
            self.sbmlprocessor.create_filter(self.work_dir + "/")
        else:
            self.sbmlprocessor.filter_filename = self.work_dir + "/" + self.sbmlprocessor.filter_filename
            self.sbmlprocessor.init()

        self.sbmlprocessor.process_sbml()
        with open(self.work_dir + 'pathways.dat', 'wb') as f:
            pickle.dump(self.sbmlprocessor.pathways, f)
        with open(self.work_dir + 'genes.dat', 'wb') as f:
            pickle.dump(self.sbmlprocessor.genes, f)
        nx.write_gpickle(self.sbmlprocessor.graph, self.work_dir + "/met.gpickle")

    def load_data(self):
        with open(self.work_dir + '/pathways.dat', 'rb') as f:
            self.sbmlprocessor.pathways = pickle.load(f)
        with open(self.work_dir + '/genes.dat', 'rb') as f:
            self.sbmlprocessor.genes = pickle.load(f)
            self.gene_reaction = defaultdict(lambda: [])
            for reaction, genes in self.sbmlprocessor.genes.items():
                for gene in genes:
                    self.gene_reaction[gene].append(reaction)
        self.sbmlprocessor.graph = nx.read_gpickle(self.work_dir + "/met.gpickle")

    def annotate_gene(self, reaction, gene_name):
        try:
            proteins = Protein.objects(organism=self.organism, gene=gene_name)
            if not len(proteins):
                proteins = Protein.objects(organism=self.organism, alias__iexact=gene_name)
            if len(proteins):
                for protein in proteins:
                    if not hasattr(protein, "reactions"):
                        protein.reactions = []  # self.add_pathways_info_to_gene(reaction, protein)
                    self.add_reaction_info_to_gene(reaction, protein)
                    self.add_properties_to_gene(protein, gene_name)
                    protein.save()
            else:

                _log.warn("%s not found" % gene_name)
                self.unmapped_genes.append(gene_name)

        except DoesNotExist:
            self.unmapped_genes.append(gene_name)
            _log.warn("%s not found" % gene_name)
        except MultipleObjectsReturned:
            self.unmapped_genes.append(gene_name)
            _log.warn(gene_name)

    def annotate_genes(self):
        genes = list(self.sbmlprocessor.graph.nodes(data=False))
        for reaction in tqdm(genes):

            if reaction in self.sbmlprocessor.genes:
                genes = self.sbmlprocessor.genes[reaction]
                if not genes:
                    continue
                for gene_name in genes:
                    if gene_name != "GENE_ASSOCIATION:":
                        self.annotate_gene(reaction, gene_name.replace("(", "").replace(")", "").replace("</p>", ""))

    def init(self):
        self.chokepoints()
        ccs = list(connected_components(self.sbmlprocessor.graph.to_undirected()))
        cp = sorted(ccs, key=lambda x: len(x))[-1]
        centrality_json = self.work_dir + "/centrality.json"
        if os.path.exists(centrality_json):
            self.centrality = json.load(open(centrality_json))
        else:
            self.centrality = {x: 0 for x in self.sbmlprocessor.graph.nodes()}
            cc = sorted(nx.connected_components(self.sbmlprocessor.graph.to_undirected()), key=lambda x: len(x))[-1]
            g2 = self.sbmlprocessor.graph.to_undirected().subgraph(cc)
            self.centrality = {x: (y if x in cp else 0) for x, y in
                               betweenness_centrality(g2).items()}
            json.dump(self.centrality, open(centrality_json, "w"))

    def annotate(self):
        if not os.path.exists(self.work_dir + "/pathways.dat"):
            self.build_network()
        else:
            _log.debug("network already calculated, loading files: " + self.work_dir + "pathways.dat")
        self.load_data()
        _log.debug("load finished")
        self.init()

        self.annotate_genes()
        self._pathways_properties()

    def pathways_items(self):
        if self.sbmlprocessor.pathways:
            return self.sbmlprocessor.pathways.items()
        else:
            return [(x["_id"], x["pathways"])
                    for x in
                    self.db.proteins.aggregate([{"$match": {"organism": self.organism, "reactions": {"$exists": 1}}},
                                                {"$project": {"reactions": 1}},
                                                {"$unwind": "$reactions"},
                                                {"$unwind": "$reactions.pathways"},
                                                {"$group": {"_id": "$reactions.name",
                                                            "pathways": {"$addToSet": "$reactions.pathways"}}}]
                                               )]

    def genes_dict(self):
        if self.sbmlprocessor.genes:
            return self.sbmlprocessor.genes
        else:
            return {x["_id"]: x["genes"] for x in
                    self.db.proteins.aggregate([{"$match": {"organism": self.organism, "reactions": {"$exists": 1}}},
                                                {"$project": {"reactions": 1, "gene": {"$arrayElemAt": ['$gene', 0]}}},
                                                {"$unwind": "$reactions"},
                                                {"$group": {"_id": "$reactions.name", "genes": {"$addToSet": "$gene"}}}
                                                ])}

    def centrality_dict(self):
        if hasattr(self, "centrality"):
            return self.centrality
        else:
            return {x["_id"]: x["c"] for x in
                    self.db.proteins.aggregate([{"$match": {"organism": self.organism, "reactions": {"$exists": 1}}},
                                                {"$project": {"reactions": 1, "c": "$search.centrality"}},
                                                {"$unwind": "$reactions"},
                                                {"$group": {"_id": "$reactions.name", "c": {"$avg": "$c"}}}]
                                               )}

    def chokepoins_dict(self):
        if hasattr(self, "_chokepoints"):
            return self._chokepoints
        else:
            return {x["_id"]: x["c"] for x in
                    self.db.proteins.aggregate([{"$match": {"organism": self.organism, "reactions": {"$exists": 1}}},
                                                {"$project": {"reactions": 1, "c": "$search.chokepoint"}},
                                                {"$unwind": "$reactions"},
                                                {"$group": {"_id": "$reactions.name", "c": {"$first": "$c"}}}]
                                               ) if x["c"]}

    def _pathways_properties(self):
        self.pathways = []
        pws_dict = {}

        genes_dict = self.genes_dict()
        chokepoints = self.chokepoins_dict()

        for _, pws in self.pathways_items():
            for pw in pws:
                pws_dict[pw] = defaultdict(list)

        for pw in pws_dict:
            pws_dict[pw]["genes"] = 0
            for r, pws in self.pathways_items():
                if pw in pws:
                    pws_dict[pw]["reactions"].append(r)
                    pws_dict[pw]["genes"] += len(genes_dict[r])
                    if genes_dict[r]:
                        pws_dict[pw]["reactions_with_gene"].append(r)

        centrality = defaultdict(lambda: 0, self.centrality_dict())
        max_centrality = max(centrality.values())

        max_reaction_pw_size = max([len(set(pws_dict[pw]["reactions"])) for pw in pws_dict])

        for pw in tqdm(pws_dict):
            reactions_with_gene = set(pws_dict[pw]["reactions_with_gene"])
            pws_dict[pw]["reactions"] = set(pws_dict[pw]["reactions"])
            pws_dict[pw]["reactions_with_gene"] = len(reactions_with_gene)
            pws_dict[pw]["chokepoints"] = len([r for r in pws_dict[pw]["reactions"] if r in chokepoints])
            try:
                pws_dict[pw]["max_centrality"] = max(
                    [centrality[r] / max_centrality for r in pws_dict[pw]["reactions"]])
            except:
                pws_dict[pw]["max_centrality"] = 0

            pws_dict[pw]["reactions"] = len(pws_dict[pw]["reactions"])
            pws_dict[pw]["norm_chokepoint"] = 1.0 * pws_dict[pw]["chokepoints"] / pws_dict[pw]["reactions"]
            pws_dict[pw]["norm_reactions"] = 1.0 * pws_dict[pw]["reactions"] / max_reaction_pw_size
            pws_dict[pw]["completeness"] = 1.0 * pws_dict[pw]["reactions_with_gene"] / pws_dict[pw]["reactions"]

            pws_dict[pw]["druggable"] = "No"
            for reaction_name in reactions_with_gene:
                druggable_reaction = self.db.proteins.count({"organism": self.organism,
                                                             "reactions.name": reaction_name,
                                                             "search.druggability": {"$gte": 0.5}})
                if druggable_reaction:
                    pws_dict[pw]["druggable"] = "Yes"
                    break

            if reactions_with_gene:
                pws_dict[pw]["human_pw_offtarget"] = 0

                for reaction_name in reactions_with_gene:
                    druggable_reaction = self.db.proteins.count({"organism": self.organism,
                                                                 "reactions.name": reaction_name,
                                                                 "search.human_offtarget": {"$lte": 0.6}})
                    if druggable_reaction:
                        pws_dict[pw]["human_pw_offtarget"] += 1

                pws_dict[pw]["human_pw_offtarget"] = (len(reactions_with_gene) - pws_dict[pw][
                    "human_pw_offtarget"]) / len(reactions_with_gene)

            ont = self.db.ontologies.find_one({"term": pw.lower()})
            if ont:
                name = ont["name"]
            else:
                name = pw

            pw_obj = PathwaySumary(term=pw, name=name, count=pws_dict[pw]["genes"], properties=pws_dict[pw])

            self.pathways.append(pw_obj)


if __name__ == "__main__":
    from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
    from SNDG.BioMongo.Process.Importer import _common_annotations, _protein_iter, import_kegg_annotation, \
        index_seq_collection, build_statistics, load_pathways

    mdb = BioMongoDB("tdr",port=27018)
    # ps = PathwaysAnnotator(mdb.db, "SaureusN315", "/data/organismos/SaureusN315/pathways/")
    # ps.sbml("Red_Staphylo_Curada_rs.sbml")
    # ps.species_filter("allfilters_con_c.dat")
    # ps.extract_genes_from_notes(lambda notes: gene_name_regexp.findall(notes))
    # ps.annotate()
    # index_seq_collection(mdb.db, "SaureusN315", pathways=True, go=True, keywords=True, ec=True, organism_idx=True,
    #                      structure=False)
    build_statistics(mdb.db, "SaureusN315")
