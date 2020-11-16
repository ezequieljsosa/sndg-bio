from Bio.KEGG import REST
from tqdm import tqdm
from time import sleep

"""
http://www.kegg.jp/kegg/docs/weblink.html
http://www.kegg.jp/kegg/rest/keggapi.html

http://rest.kegg.jp/list/pathway
http://rest.kegg.jp/get/hsa00650/kgml
http://rest.kegg.jp/get/ko00010/kgml --> generico
http://www.kegg.jp/kegg-bin/show_pathway?@ko00010/reference%3dwhite/default%3d%23bfffbf/K00844/K01810/K00850/K00895/K03841/K01623/K01803/K00134/K00927/K01834/K15633/K15634/K01689/K00873/K00161/K00162/K00627/K00382/K00016/K01568/K00121/K18857/K00001/K00002/K00128/K14085/K01895/K01785/K01835/K01792/K03103/K00131/K01610

http://www.kegg.jp/kegg-bin/download_htext?htext=br08901.keg&format=json&filedir=

"""

# for db in ["pathway","ko","cpd","brite"]:
#     with open("/data/databases/kegg/" + db + ".txt","w") as h:
#         print db
#         data = REST.kegg_list(db).read()
#         print len(data)
#         h.write(data)
# wget http://www.kegg.jp/kegg-bin/download_htext?htext=br08901.keg&format=json&filedir=
# wget http://www.kegg.jp/kegg-bin/download_htext?htext=br08001.keg&format=json&filedir=

# L = list(open("/data/databases/kegg/pathway.txt"))
#
# for pathway in tqdm(L):
#     pw = "ko" + pathway.split()[0].split(":map")[1]
#     with open("/data/databases/kegg/ko/" + pw + ".kgml", "w") as h:
#         print pw
#         try:
#             data = REST.kegg_get(pw, option="kgml").read()
#             h.write(data)
#             sleep(1)
#         except:
#             pass

import glob
import os
from collections import defaultdict
from Bio.KEGG.KGML.KGML_parser import read


# for x in glob.glob("/data/databases/kegg/ko/*.kgml"):
#     if os.path.getsize(x) > 100:
#         print x
#         pepe = read(open("/data/databases/kegg/ko/ko00010.kgml"))

# print

# print REST.kegg_get("ko:K00001","json")

class Kegg(object):
    """

    """

    @staticmethod
    def update_files(base_dir="/data/databases/kegg/"):
        for db in ["pathway","ko","cpd","brite"]:
            with open(base_dir + db + ".txt","w") as h:
                data = REST.kegg_list(db).read()

                h.write(data)
        # wget http://www.kegg.jp/kegg-bin/download_htext?htext=br08901.keg&format=json&filedir=
        # wget http://www.kegg.jp/kegg-bin/download_htext?htext=br08001.keg&format=json&filedir=

        L = list(open(base_dir + "pathway.txt"))

        for pathway in tqdm(L):
            pw = "ko" + pathway.split()[0].split(":map")[1]
            kgmlpath = base_dir + "ko/" + pw + ".kgml"
            if not os.path.exists(kgmlpath):
                with open(kgmlpath, "w") as h:

                    try:
                        data = REST.kegg_get(pw, option="kgml").read()
                        h.write(data)
                        sleep(1)
                    except:
                        pass

    def __init__(self, ko_list_path="/data/databases/kegg/ko.txt",
                 brittle_pw_path="/data/databases/kegg/ko/",
                 kgmls_dir="/data/databases/kegg/ko/"):
        self.ko_list_path = ko_list_path
        self.brittle_pw_path = brittle_pw_path
        self.kgmls_dir = kgmls_dir

        self.ko_dict = defaultdict(lambda :{})
        self.pw_dict = {}
        self.ko_react = {}


    def init(self):

        with open(self.ko_list_path) as h:
            for koline in h.readlines():
                ko, desc = koline.replace("#160;","").split("\t")
                gene = desc.split(";")[0]
                desc2 = ";".join( desc.split(";")[1:])
                desc = desc2.strip().split("[EC:")[0]
                ecs = None
                if len(desc2.split("[EC:")) > 1:
                    ecs = desc2.split("[EC:")[1].replace("]","")
                    ecs = ["ec:" + x for x in ecs.split()]
                self.ko_dict[ko]["genes"] = gene.split()
                self.ko_dict[ko]["ecs"] = ecs
                self.ko_dict[ko]["desc"] = desc
                self.ko_dict[ko]["pathways"] = []

        for kgml_path in glob.glob(self.kgmls_dir + "/*.kgml"):
            if os.path.getsize(kgml_path) > 100:
                kgml = read(open(kgml_path))
                pw_name = "kegg_pw:" + kgml.name.split(":ko")[1]

                ko_react = defaultdict(lambda: [])

                for  ortholog in kgml.orthologs:

                    for ko in ortholog.name.split():

                        if ":" in ko and "pathways" in self.ko_dict[ko]:
                            self.ko_dict[ko]["pathways"].append(pw_name)
                            if ortholog.reaction:
                                ko_react[ko.split(":")[1]].append(ortholog.reaction)

                self.pw_dict[pw_name] = {
                    "name":pw_name,
                    "title": kgml.title,
                    "reactions": dict(ko_react)
                }
        self.ko_dict = dict(self.ko_dict)

    def read_annotation(self, query_ko_path):
        self.ko_gene = defaultdict(lambda: [])
        with open(query_ko_path) as h:
            for x in h:
                x = x.strip().split("\t")
                if len(x) > 1:
                    self.ko_gene[x[1]].append(x[0])




if __name__ == '__main__':
    kegg = Kegg(ko_list_path="/data/databases/kegg/ko.txt",
                brittle_pw_path="/data/databases/kegg/ko/",
                kgmls_dir="/data/databases/kegg/ko/")
    kegg.init()
    ilex_data = "/data/organismos/ILEX_PARA_TRANSCRIPT/annotation/query.ko"
    kegg.read_annotation(ilex_data)

