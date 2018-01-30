"""

"""
import os
import glob
import logging
from tqdm import tqdm
import json

import pandas as pd
import Bio.SearchIO as bpsio

_log = logging.getLogger(__name__)


class Modelome(object):

    columns = [ "model",
      "prot",
      "template",
      "identity",
      "zqmean",
      "ga341",
      "dope",
      "qstart",
      "qend",
      "hstart",
      "hend",
      "hresstart",
      "qaln",
      "haln",
      "evalue"]

    def load_result_from_dir(self, base_path):
        """
        Creates a csv with all the built models in the directory
        :return:
        """

        directories = [o for o in os.listdir(base_path)
                       if os.path.isdir(o)]

        no_hit = 0
        no_model = 0
        best = {}
        with tqdm(directories) as pbar:
            for s in pbar:
                pbar.set_description("Processing %s" % s.id)
                if os.path.getsize("./" + s.id + "/profile_search.xml"):
                    quals = []
                    for m in glob.glob("./" + s.id + "/" + s.id + "/*"):
                        for mi in glob.glob(m + "/*.B9999000?.pdb.json"):
                            with open(mi) as h:
                                q = json.load(h)["z_score"]
                                quals.append([mi.replace(".json", ""), q])
                    if quals:
                        x = sorted(quals, key=lambda x: x[1])
                        best[s.id] = sorted(quals, key=lambda x: abs(x[1]))[0]
                    else:
                        no_hit += 1
                else:
                    no_model += 1
        self.result = {"no_hit": no_hit, "no_model": no_model, "models": best}

    def build_structures_table(self):
        self.df = pd.DataFrame()
        for prot, (model_path, zqmean) in tqdm(self.result["models"].items()):
            with open(model_path) as h2:
                lines = h2.readlines()
                data = [x.replace("REMARK   6 ", "") for x in lines if x.startswith("REMARK")]
            data = {x.split(":")[0]: x.split(":")[1].strip() for x in data}
            aln = bpsio.read("./" + prot + "/profile_search.xml", "blast-xml")
            template = model_path.split(prot)[-2][1:-1]
            hsp = [hit[0] for hit in aln if hit.id == template][0]

            struct = {
                "model": model_path.split(os.sep)[-2],
                "prot": prot,
                "template": template,
                "identity": data["MODELLER BEST TEMPLATE % SEQ ID"],
                "zqmean": abs(zqmean),
                "ga341": float(data["GA341 score"]),
                "dope": float(data["DOPE score"]),
                "qstart": int(hsp.query_start),
                "qend": int(hsp.query_end),
                "hstart": int(hsp.hit_start),
                "hend": int(hsp.hit_end),
                "hresstart": int(data["TEMPLATE"].split()[1].split(":")[0]),
                "qaln": str(hsp.aln[0].seq),
                "haln": str(hsp.aln[1].seq),
                "evalue": hsp.evalue
            }
            self.df = self.df.append(struct,ignore_index=True)
