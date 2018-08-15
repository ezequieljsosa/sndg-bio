"""

"""
import glob
import json
import logging
import os
import shutil
from collections import defaultdict

import pandas as pd

from tqdm import tqdm

import Bio.SearchIO as bpsio
from SNDG import Struct
from SNDG import mkdir
from SNDG.Structure.ChainSplitter import ChainSplitter
from SNDG.Structure.Modeller import Modeller


_log = logging.getLogger(__name__)

from modeller.automodel import assess, refine  # @UnresolvedImport

REFINEMENT = refine.very_slow
ASSESMENTS = [assess.GA341, assess.DOPE]  # @UndefinedVariable
MODELS_TO_GENERATE = 1


def default_hsp_fn(base_path, prot, template):
    aln = bpsio.read(base_path + "/" + prot + "/profile_search.xml", "blast-xml")
    hsp = [hit[0] for hit in aln if hit.id == template][0]
    return hsp


class Modelome(object):
    columns = ["model",
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

    @staticmethod
    def load_result_from_dir(base_path,
                             model_fn=lambda base_path, protein_name: glob.glob(
                                 base_path + "/" + protein_name + "/" + protein_name + "/*")):
        """
        Creates a csv with all the built models in the directory
        :return:
        """

        directories = [o for o in os.listdir(base_path)
                       if os.path.isdir(base_path + "/" + o) and not o.startswith(".")]

        no_hit = 0
        no_model = 0
        best = {}
        with tqdm(directories) as pbar:
            for s in pbar:
                pbar.set_description("Processing %s" % s)
                models = model_fn(base_path,s)
                if models:
                    quals = []
                    for m in models:
                        for mi in glob.glob(m + "/*.B9999000?.pdb.json"):
                            with open(mi) as h:
                                q = json.load(h)["z_score"]
                                quals.append([mi.replace(".json", ""), q])
                    if quals:
                        best[s] = sorted(quals, key=lambda xx: abs(xx[1]))[0]
                    else:
                        no_hit += 1
                else:
                    no_model += 1
        return {"no_hit": no_hit, "no_model": no_model, "models": best}

    @staticmethod
    def build_structures_table(base_path, model_zqmean_tuples, hsp_fn=default_hsp_fn):
        df = pd.DataFrame()
        for prot, (model_path, zqmean) in tqdm(model_zqmean_tuples):
            with open(model_path) as h2:
                lines = h2.readlines()
                data = [x.replace("REMARK   6 ", "") for x in lines if x.startswith("REMARK")]
            data = {x.split(":")[0]: x.split(":")[1].strip() for x in data}
            template = model_path.split(prot)[-2][1:-1]
            hsp = hsp_fn(base_path, prot,template)
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
                "evalue": hsp.evalue,
                "path": model_path
            }
            df = df.append(struct, ignore_index=True)
        return df

    @staticmethod
    def prepare_dir(directory, df, mfilter=lambda df: df[(df.zqmean >= -2) & (df.zqmean <= 2)],
                    csv_name="models.csv"):
        df = mfilter(df)
        mkdir(directory)
        df.to_csv(directory + "/" + csv_name, index=False, columns=Modelome.columns)
        for _, r in df.iterrows():
            shutil.copy(r.path, directory + "/" + r.model + ".pdb")

    @staticmethod
    def model_hsps(seq_id, work_dir, hsps, refinement=REFINEMENT, models_to_generate=MODELS_TO_GENERATE,
                   assessments=ASSESMENTS, entries={}, pdb_divided="/data/databases/pdb/divided/",
                   tmp_dir=None,max_models=3):
        result  = {"models":defaultdict(lambda: {})}

        alns = []
        for hsp in hsps:
            aln = Struct(
                aln_query=Struct(name=seq_id, seq=str(hsp.aln[0].seq), start=hsp.query_start,
                                 end=hsp.query_end),
                aln_hit=Struct(name=hsp.hit.id, seq=str(hsp.aln[1].seq), start=hsp.hit_start, end=hsp.hit_end))
            alns.append(aln)

        modeler = Modeller(work_dir, tmp_dir)

        modeler._refinement = refinement
        modeler.model_count = models_to_generate
        modeler._assess_methods = assessments
        modeler.parallel_jobs = 1

        def pdb_fn(x):
            return x.aln_hit.name.split("_")[0]

        alns = sorted(alns, key=lambda x: entries[pdb_fn(x)] if pdb_fn(x) in entries else 20)
        result["alns"] = alns
        for aligment in alns[0:max_models]:
            # pdb,aligment = pdb_alignment
            pdb, chain, _, _ = aligment.aln_hit.name.split("_")
            if not os.path.exists(modeler.pdb_path(seq_id + "_" + aligment.aln_hit.name, seq_id)):
                base_model_path = pdb_divided + pdb[1:3] + "/pdb" + pdb + ".ent"
                ChainSplitter(tmp_dir).make_pdb(base_model_path, pdb, chain, overwrite=True)
                models = modeler.create_model(seq_id + "_" + aligment.aln_hit.name, aligment)
            else:
                models = [modeler.pdb_path(seq_id + "_" + aligment.aln_hit.name, seq_id, idx) for idx in
                          range(1, models_to_generate + 1)]
            result["models"][aligment.aln_hit.name] = models
        return result



if __name__ == '__main__':
    from SNDG import execute
    from SNDG.Structure.QuickModelome import QuickModelome


    result = Modelome.load_result_from_dir("/data/organismos/GCF_001624625.1/structures/results/")
    df = Modelome.build_structures_table("/data/organismos/GCF_001624625.1/structures/results/", result["models"].items())
    Modelome.prepare_dir("/data/organismos/GCF_001624625.1/structures/final", df,mfilter=lambda x:x)

    # m = Modelome()
    # qm = QuickModelome()
    # qm.load_hsp_dict("/data/pext/good_hits.xml")
    # result = m.load_result_from_dir("/data/pext", model_fn=QuickModelome.protein_models)
    # df = m.build_structures_table("/data/pext", result["models"].items(), hsp_fn=qm.aln_hsp)
    # m.prepare_dir("/data/pext_final", df, csv_name="models2.csv")
