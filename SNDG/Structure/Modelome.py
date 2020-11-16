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

from SNDG.Structure.Modeller import Modeller, ModellerOverflowError

_log = logging.getLogger(__name__)

from modeller.automodel import assess, refine  # @UnresolvedImport
from modeller import SequenceMismatchError
from _modeller import ModellerError

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
                models = model_fn(base_path, s)
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
            hsp = hsp_fn(base_path, prot, template)
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
                   assessments=ASSESMENTS, entries={},
                   tmp_dir=None, max_models=3):
        result = {"models": defaultdict(lambda: {}), "errors": []}

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
            if ";" in aligment.aln_hit.name:
                pdb = aligment.aln_hit.name[-5:-1]
                chain = aligment.aln_hit.name[-1]
                modeller_pdb_id = f'{seq_id}_{pdb}_{chain}'
                aligment.aln_hit.name = f'{pdb}_{chain}'
            else:
                pdb, chain, _, _ = aligment.aln_hit.name.split("_")
                modeller_pdb_id = f'{seq_id}_{pdb}_{chain}'

            if not os.path.exists(modeler.pdb_path(modeller_pdb_id, seq_id)):
                try:
                    models = modeler.create_model(modeller_pdb_id, aligment)
                except SequenceMismatchError as ex:
                    result["errors"].append(f"modeller_pdb_id {str(ex)}\n")
                    continue
                except ModellerError as ex:
                    result["errors"].append(f"modeller_pdb_id {str(ex)}\n")
                    continue

            else:
                models = [modeler.pdb_path(modeller_pdb_id, seq_id, idx) for idx in
                          range(1, models_to_generate + 1)]
            result["models"][aligment.aln_hit.name] = models
        result["models"] = dict(result["models"])
        return result

    @staticmethod
    def process_file(params):
        aln_file = params["aln_file"]
        templates2use = params["templates2use"]
        tmp_dir = params["tmp_dir"]
        output_dir = params["output_dir"]

        aln_file = aln_file.strip()
        if os.path.getsize(aln_file) < 100:
            return [{"errors": f'\n{aln_file} empty file'}]
        hsps = []
        try:
            hsps = [hsp
                    for query_result in bpsio.parse(aln_file.strip(), "blast-xml")
                    for hit in query_result
                    for hsp in hit]
        except ValueError:
            sys.stderr.write(f"error reading alignments in {aln_file}")

        hsps = hsps[:templates2use]
        if hsps:
            seq_id = hsps[0].query.id
            # pdb_chains = [x.split("_") for x in set([hsp.hit.id[3:7] + "_" + hsp.hit.id[-1] for hsp in hsps])]
            pdb_chains = [[hsp.hit.id[3:7], hsp.hit.id[-1]] for hsp in hsps]
            updated = False
            for pdb, _ in pdb_chains:
                pdb_utils.update_pdb(pdb)
                updated =  os.path.exists(pdb_utils.pdb_path(pdb))
                break

            if not updated:
                sys.stderr.write(f'{pdb} could not be updated...\n')
                return
            pdb_utils.extract_chains(pdb_chains, tmp_dir)
            models_results = []
            for hsp in hsps:
                try:
                    models_result = Modelome.model_hsps(seq_id, os.path.abspath(output_dir), [hsp], refinement=REFINEMENT,
                                                models_to_generate=MODELS_TO_GENERATE,
                                                assessments=ASSESMENTS, entries={},
                                                tmp_dir=tmp_dir, max_models=1)
                except ModellerOverflowError as e:
                    sys.stderr.write(f"error processing {seq_id}: {str(e)}")
                    continue
                models_results.append(models_result)
            return models_results

        else:
            return [{"errors": f'\nno aligments for {aln_file}\n'}]


if __name__ == "__main__":
    import argparse
    import sys
    import Bio.SeqIO as bpio
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from SNDG import mkdir
    from SNDG import execute
    from tempfile import mkdtemp
    from SNDG.Structure.PDBs import PDBs
    from tqdm.contrib.concurrent import process_map

    parser = argparse.ArgumentParser(description='Profile utils')

    parser.add_argument('alns', type=argparse.FileType('r'),
                        default=sys.stdin, help='file with list of file alignments to process.')
    """Formats: .xml(blast outfmt 5), fasta or table"""
    parser.add_argument("-o", '--output_dir', default="./", help='generated models dir')
    parser.add_argument("-d", '--pdbs_dir', default="/data/databases/pdb/", help='pdb files directory')
    parser.add_argument('--tmp_dir', default=mkdtemp(), help='temporal directory')
    parser.add_argument('--cpus', default=1, type=int, help='cpu cores to use')
    # parser.add_argument('--max_alns', default=3, type=int, help='max different templates to use')
    parser.add_argument('-t', "--templates_to_use", default=3, type=int, help='max amount of templates to use.')

    args = parser.parse_args()
    pdbs_dir = args.pdbs_dir + ("/" if args.pdbs_dir[-1] != "/" else "")
    mkdir (f'{pdbs_dir}/divided')
    pdb_utils = PDBs(pdbs_dir)
    # pbar = tqdm(args.alns)
    sys.stderr.write(str(args))
    sys.stderr.write(f'reading alignment file\n')
    alns = [{"aln_file": x, "templates2use": args.templates_to_use,
             "output_dir": args.output_dir, "tmp_dir": args.tmp_dir} for x in args.alns]
    mkdir(args.output_dir)
    assert os.path.exists(args.output_dir),f'"{args.output_dir}" could not be created'

    sys.stderr.write(f'processing alignment files\n')
    for result in process_map(Modelome.process_file,
                              alns,  max_workers=args.cpus, chunksize=1):
        #,description="alignment files to process"
        print(result)

        # pbar.set_description(f"processing: {aln_file}")

    # result = Modelome.load_result_from_dir("/data/organismos/GCF_001624625.1/structures/results/")
    # df = Modelome.build_structures_table("/data/organismos/GCF_001624625.1/structures/results/", result["models"].items())
    # Modelome.prepare_dir("/data/organismos/GCF_001624625.1/structures/final", df,mfilter=lambda x:x)

    # m = Modelome()
    # qm = QuickModelome()
    # qm.load_hsp_dict("/data/pext/good_hits.xml")
    # result = m.load_result_from_dir("/data/pext", model_fn=QuickModelome.protein_models)
    # df = m.build_structures_table("/data/pext", result["models"].items(), hsp_fn=qm.aln_hsp)
    # m.prepare_dir("/data/pext_final", df, csv_name="models2.csv")
