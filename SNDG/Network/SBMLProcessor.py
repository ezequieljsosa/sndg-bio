'''
Created on Jan 9, 2015

@author: eze

To install
sudo apt install libsbml5-dev
pip install python-libsbml

'''

import os
import re
import logging
from tqdm import tqdm
import networkx as nx
from libsbml import readSBML

_log = logging.getLogger("SBMLProcessor")


class SBMLProcessor(object):
    '''
    Process the output of a PathwayTools generated sbml file
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.ignore_filter = False;
        self.ignore_reversibility = False;
        self.sbml_file_path = None
        self.filter_filename = "allfilters_con_c.dat";
        import re
        reg = r'\((.*?)\)'
        # self.fn_extract_genes = lambda x: [y.strip().replace("(", "").replace(")", "").replace("</p>", "") for y in
        #                                    x.split() if y.strip() and y.strip() != 'GENE_ASSOCIATION:' and y != "or"]
        self.fn_extract_genes = lambda x: re.findall(reg, x)
        self.graph = nx.DiGraph()
        self.pathways = {}
        self.genes = {}
        self.filter_list = []

    def decode_sbml(self, text):
        result = text
        replacements = [("__45__", "-"), ("__46__", "."), ("__43__", "+"), ("_47__", "/"),
                        ("_1", "1"), ("_2", "2"), ("_3", "3")
            , ("_4", "4"), ("_5", "5"), ("_6", "6"), ("_7", "7"), ("_8", "8"), ("_9", "9")]
        for original, replacement in replacements:
            result = result.replace(original, replacement)
        result = re.sub(r"_c$", "", result)
        return result

    def complete_reaction_data(self, reaction):
        rxnId = self.decode_sbml(reaction.getName())
        notas = reaction.getNotesString()
        pathways, genes = self.extract_pathways_and_genes(rxnId, notas)
        self.genes[rxnId] = genes
        self.pathways[rxnId] = pathways

    def extract_pathways_and_genes(self, rxnId, notas):

        self.genes[rxnId] = None

        genes, pathways = ([], [])
        if (notas.find("GENE_ASSOCIATION") != -1 or
                notas.find("SUBSYSTEM") != -1):  # possui Rv ou pathway preenchido
            listaNotas = notas.split("<p>")
            for l in listaNotas:
                if l.startswith("GENE_ASSOCIATION"):
                    locus = self.fn_extract_genes(l)  # retorna todos os matches com a ER
                    genes = locus
                if l.startswith("SUBSYSTEM"):
                    pathways = l[l.find(":") + 1:l.rfind("<")].strip().split(",")
                    pathways = [self.decode_sbml(k.strip()) for k in pathways]
        return pathways, genes

    def init(self):
        assert self.sbml_file_path, "SBML file not set"
        self.document = readSBML(self.sbml_file_path);
        if os.path.exists(self.filter_filename):
            with open(self.filter_filename, "r") as filter_file:
                self.filter_list = [x.split("\t")[0].strip() for x in filter_file.read().lower().splitlines()]
        else:
            self.filter_list = []

        # errors = self.document.getNumErrors();
        # if (errors > 0):
        #    _log.error(str(errors))
        #    raise ("Encountered the following SBML errors:" + str(errors))

        self.model = self.document.getModel();

    def create_filter(self, filter_path):
        metabolites_count = {}
        for reaction in self.model.getListOfReactions():
            for metabolite in set([x.species for x in reaction.getListOfProducts()] + [x.species for x in
                                                                                       reaction.getListOfReactants()]):
                if metabolite in metabolites_count:
                    metabolites_count[metabolite] += 1
                else:
                    metabolites_count[metabolite] = 1
        limit = 20
        with open(filter_path + self.filter_filename, "w") as h:
            for met, count in reversed(sorted(metabolites_count.items(), key=lambda x: x[1])):
                if count >= limit:
                    self.filter_list.append(self.decode_sbml(met.lower()))
                    h.write(self.decode_sbml(met.lower()) + "\t" + str(count) + "\n")

    def add_reaction_to_graph(self, reaction):
        if not self.graph.has_node(self.decode_sbml(reaction.getName())):
            products = [self.decode_sbml(x.getSpecies()) for x in reaction.getListOfProducts() if
                        self.decode_sbml(x.getSpecies()) not in self.filter_list]
            rectants = [self.decode_sbml(x.getSpecies()) for x in reaction.getListOfReactants() if
                        self.decode_sbml(x.getSpecies()) not in self.filter_list]
            self.graph.add_node(self.decode_sbml(reaction.getName()), products=products, reactants=rectants)

    def process_sbml(self):

        listOfReactions = self.model.getListOfReactions();

        for reaction in tqdm(listOfReactions):
            self.complete_reaction_data(reaction)
            # get list of Products (and Reactants if the reaction is reversible)
            listOfReactantsAndProducts = reaction.getListOfProducts().clone();

            if self.ignore_reversibility or reaction.getReversible():
                listOfReactantsAndProducts.appendFrom(reaction.getListOfReactants());

            self.add_reaction_to_graph(reaction)

            for element1 in listOfReactantsAndProducts:
                elem1_species = element1.getSpecies();
                if (self.ignore_filter or (self.decode_sbml(elem1_species).lower() not in self.filter_list)):
                    listOfReactions2 = listOfReactions.clone();
                    listOfReactions2.remove(reaction.id);  # remove reaction1 to avoid checking against itself
                    for reaction2 in listOfReactions2:

                        self.add_reaction_to_graph(reaction2)

                        # get list of Reactants (and Products if the reaction is reversible)
                        listOfReactantsAndProducts2 = reaction2.getListOfReactants().clone();
                        if self.ignore_reversibility or reaction2.getReversible():
                            listOfReactantsAndProducts2.appendFrom(reaction2.getListOfProducts());

                        for element2 in listOfReactantsAndProducts2:
                            elem2_species = element2.getSpecies();
                            if self.ignore_filter or (self.decode_sbml(elem2_species).lower() not in self.filter_list):
                                if element1.species == element2.species:
                                    self.graph.add_edge(self.decode_sbml(reaction.getName()),
                                                        self.decode_sbml(reaction2.getName()),
                                                        linkedby=self.decode_sbml(element1.species))

    def toSIF(self, filename, cc=True):

        if cc:
            g = self.graph.subgraph(
                sorted([(x, len(x)) for x in nx.connected_components(self.graph.to_undirected())], key=lambda y: y[1])[0][0])
        else:
            g = self.graph

        with open(filename, "w") as filehandler:
            for reaction, reaction2 in g.edges:
                link = reaction + "\tlinkedWith\t" + reaction2 + "\n";
                filehandler.write(self.decode_sbml(link));
                # else:
                #     print  >> filehandler, self.decode_sbml(reaction.id);


"""Reactions graph
from networkx.algorithms import bipartite
document = readSBML("./pathways-sm.sbml");
model = document.getModel() 

reacts = []                               
comps = []           
edges = []                                                                                
for reaction in model.getListOfReactions():                  
    r = decode_sbml(reaction.getName())                                                    
    reacts.append(  r )                           
    for s in reaction.getListOfProducts():                
        s = decode_sbml(s.getSpecies())    
        if s.lower() not in exclude:
            comps.append(s)               
            edges.append( (r,s,) )         
    for s in reaction.getListOfReactants():
        s = decode_sbml(s.getSpecies())
        if s.lower() not in exclude:
            comps.append(s)
            edges.append( (s,r,) )
                 
g = nx.DiGraph()
g.add_nodes_from(reacts, bipartite=0)
g.add_nodes_from(comps, bipartite=1)
g.add_edges_from(edges)


production = []
consuption = []
for x in comps:
    if g.out_degree(x) == 1:
        e = list(g.out_edges(x))[0]
        consuption.append(e[1])
    if g.in_degree(x) == 1:
        e = list(g.in_edges(x))[0]
        production.append(e[0])
        
double = set(production) & set(consuption)
production = set(production) - double
consuption = set(consuption) -double
print("Double")
print(double)
print("Prod")
print(production)
print("Consu")
print(consuption)

"""

if __name__ == "__main__":
    import argparse
    import os
    import sys

    parser = argparse.ArgumentParser(description='SBML utils')
    parser.add_argument('-i', '--sbml', help="sbml", required=True)
    parser.add_argument('-f', '--filter', help="List of compounds to be filtered. One per line",
                        default="ubiquitous_compounds.txt")
    parser.add_argument('-o', '--output_dir', help="output dir", default="./")

    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    assert os.path.exists(args.output_dir), f'{args.output_dir} could not be created'

    sbml = SBMLProcessor()
    sbml.sbml_file_path = args.sbml
    sbml.filter_filename = args.filter
    if os.path.exists(f"{args.output_dir}/network.gpickle"):
        sbml.graph = nx.read_gpickle(f"{args.output_dir}/network.gpickle")
    else:
        sbml.init()
        if not os.path.exists(args.filter):
            sbml.create_filter(args.output_dir)
            print("filter fille was just created, check it and re run the command with: ")
            print(f" -i '{args.sbml}' -o '{args.output_dir}' -f '{args.filter}'")
            sys.exit(0)
        sbml.process_sbml()
        nx.write_gpickle(sbml.graph, args.output_dir + "/network.gpickle")
    sbml.toSIF(f"{args.output_dir}/network.sif",cc=False)
