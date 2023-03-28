from collections import defaultdict

import pythoncyc
import pickle
import networkx as nx
import subprocess as sp
from threading import Thread

def execute(cmd):
    sp.check_output(cmd,shell=True)
    """
    popen = sp.Popen(cmd, stdout=sp.PIPE, universal_newlines=True, shell=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise sp.CalledProcessError(return_code, cmd)
    """
class PathoLogic:
    DOCKERIMAGENAME = "ezequieljsosa/pathwaytools:2.4"
    DOCKERCONTAIERNAME = "pathwaytools"

    def __init__(self, orgdbname, pgdbs_data_dir, input_data_dir, output_data_dir, filter_count=20):

        self.orgdbname = orgdbname
        self.pgdbs_data_dir = pgdbs_data_dir
        self.input_data_dir = input_data_dir
        self.output_data_dir = output_data_dir
        self.server_thread = None
        self.filter_count = filter_count

    def run_server(self):
        print("jhgjksdfghksdjlfhg")
        execute(
            f'docker run --rm --name {PathoLogic.DOCKERCONTAIERNAME} -v {self.input_data_dir}:{self.input_data_dir} \
                    -w {self.input_data_dir} --volume {self.output_data_dir}:/opt/data/ptools-local/pgdbs \
                    -p 5008:5008 {PathoLogic.DOCKERIMAGENAME} /opt/pathway-tools/pathway-tools -python -api')
        print("aaaaaaajhgjksdfghksdjlfhg")

    def start(self):
        execute(f'docker start {PathoLogic.DOCKERCONTAIERNAME}')

    def stop(self):
        execute(f'docker start {PathoLogic.DOCKERCONTAIERNAME}')

    def run_pathologic(self):
        execute(f'docker exec -w {self.input_data_dir} {PathoLogic.DOCKERCONTAIERNAME} /opt/pathway-tools/pathway-tools \
                     -no-patch-download -no-cel-overview -no-web-cel-overview -patho {self.input_data_dir}')

    def process_pgdb(self):
        self.db = pythoncyc.select_organism(self.orgdbname)
        # pathwayts = {x.replace("|",""): self.db.get_frame_objects([x])[0] for x in self.db.all_pathways()}
        reactions = {x.replace("|", ""): self.db.get_frame_objects([x])[0] for x in self.db.all_reactions()}
        # compounds = {x["frameid"].replace("|",""): self.db.get_frame_objects([x["frameid"]])[0] for x in self.db.compounds}
        genes = {x["frameid"].replace("|", ""): self.db.reactions_of_gene(x["frameid"]) for x in self.db.genes}

        genes2 = {}

        graph = nx.DiGraph()
        with open(self.output_data_dir + "/pathways.dat", "wb") as hp, open(
                self.output_data_dir + "/genes.dat", "wb"
        ) as hg, open(self.output_data_dir + "/met.gpickle", "wb") as hr:

            # pathwayts2 = { k:v["common_name"] for k,v in pathwayts.items()}
            # if x in  ["common_name", "names", "score"]}
            pathwayts3 = {}
            links_prod = defaultdict(list)
            links_sub = defaultdict(list)

            for k, r in reactions.items():
                if "in_pathway" in r.__dict__:
                    pathwayts3[k] = [p.replace("|", "") for p in r["in_pathway"]]
                else:
                    pathwayts3[k] = []
                products = [x.replace("|", "") for x in r.right if x]
                reactans = [x.replace("|", "") for x in r.left if x]
                genes2[k] = []
                graph.add_node(k, products=products, reactants=reactans)
                for p in products:
                    links_prod[p].append(k)
                    links_sub[p].append(k)

            pickle.dump(pathwayts3, hp)

            for gene, reactions2 in genes.items():
                if reactions2:
                    for reaction in reactions2:
                        if reaction[1:-1] in genes2:
                            genes2[reaction[1:-1]].append(gene)
                        else:
                            genes2[reaction[1:-1]] = []

            for comp in (set(links_sub) & set(links_prod)):
                count = len(links_sub[comp]) + len(links_prod[comp])
                if count > self.filter_count:
                    # print( comp + ": " +  str(count))
                    continue
                for rsub in links_sub[comp]:
                    for rprod in links_prod[comp]:
                        graph.add_edge(rprod, rsub, linkedby=comp)

            pickle.dump(graph, hr)
            pickle.dump(genes2, hg)




if __name__ == "__main__":
    import argparse
    import os
    from time import sleep

    parser = argparse.ArgumentParser(description='Pathologic')
    # orgdbname,pgdbs_data_dir,input_data_dir,output_data_dir,filter_count
    parser.add_argument('orgdbname', help="name if the pgdb in pathwaytools")
    parser.add_argument('pgdbs_data_dir', help="absolute path to mount PGDBs directory")
    parser.add_argument('input_data_dir',
                        help="directory with Pathologic input (organism, genetic elements and genebank)")
    parser.add_argument('output_data_dir', help="directory where results are stores")
    parser.add_argument('--filter_count', type=int, help="max occurrences for compounds")
    parser.add_argument('--wait_server_up', type=int, default=10,
                        help="time to wait before using pathwaytools servere")

    args = parser.parse_args()

    assert os.path.exists(args.pgdbs_data_dir), f'{args.pgdbs_data_dir} does no exist'
    assert os.path.exists(args.input_data_dir), f'{args.input_data_dir} does no exist'
    if not os.path.exists(args.output_data_dir):
        os.makedirs(args.output_data_dir)
    assert os.path.exists(args.output_data_dir), f'{args.output_data_dir} cant be created'

    pl = PathoLogic(args.orgdbname, args.pgdbs_data_dir, args.input_data_dir, args.output_data_dir, args.filter_count)

    pl.server_thread = Thread(target=pl.run_server)
    pl.server_thread.start()
    try:
        sleep(args.wait_server_up)
        pl.run_pathologic()
        pl.stop()
        pl.server_thread.join()
        pl.server_thread = Thread(target=pl.run_server)
        pl.server_thread.start()
        sleep(args.wait_server_up)
        pl.process_pgdb()
        pl.stop()
        pl.server_thread.join()
    finally:
        pl.stop()
