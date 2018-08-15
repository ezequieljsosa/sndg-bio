"""

"""

import glob
import logging
from collections import defaultdict

from tqdm import tqdm

import Bio.SearchIO as bpsio
from SNDG.Structure.Modelome import Modelome

import argparse
import multiprocessing


_log = logging.getLogger(__name__)


class QuickModelome(object):
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

    def __init__(self):
        self.hsp_dict = defaultdict(lambda: {})

    def load_hsp_dict(self, xml_blast_result):

        for query in bpsio.parse(xml_blast_result, "blast-xml"):
            for hit in query:
                if list(hit):
                    hsp = list(hit)[0]
                    self.hsp_dict[query.id][hsp.hit.id] = hsp

    @staticmethod
    def protein_models(work_dir, protein_name):
        return glob.glob(work_dir + "/" + protein_name + "/" + "/*")

    def aln_hsp(self, base_path, prot, template):
        return self.hsp_dict[prot][template]

    def quick_structurome(self, xml_blast_result, data_dir, entries, tmp_dir="/tmp/chain_PDBs",
                          pdb_divided="/data/databases/pdb/divided/", max_models=3):

        good_model = defaultdict(lambda: [])

        def identity(hsp):
            return 1.0 * hsp.ident_num / hsp.aln_span

        _log.info("searching good templates")
        for query in tqdm(bpsio.parse(xml_blast_result, "blast-xml")):
            for hit in query:
                if list(hit):
                    hsp = list(hit)[0]
                    if 0.6 <= identity(hsp) < 0.95:
                        good_model[hsp.query.id].append(hsp)



        tuplas = good_model.items()

        _log.info("creating models")
        with tqdm(tuplas) as pbar:
            for seq, hsps in pbar:
                try:
                    from SNDG.Structure.Modelome import Modelome
                    Modelome.model_hsps(seq, data_dir, hsps, entries=entries, tmp_dir=tmp_dir,
                                        pdb_divided=pdb_divided, max_models=max_models)
                except Exception as ex:
                    _log.exception(ex)


if __name__ == '__main__':
    from SNDG import init_log,mkdir
    init_log()
    import os
    parser = argparse.ArgumentParser()

    parser.add_argument("-cpus", "--cpus", default=multiprocessing.cpu_count())
    parser.add_argument("--xml_blast_result", required=True)
    parser.add_argument("--data_dir")
    parser.add_argument("--entries", default='/data/pdb/entries.idx')
    parser.add_argument("--pdb_divided", default="/data/pdb/divided/")
    parser.add_argument("--max_models",  type=int,  default=3)
    parser.add_argument("--tmp_dir",   default="/tmp/chain_PDBs")

    args = parser.parse_args()
    assert os.path.exists(args.tmp_dir),args.tmp_dir
    assert os.path.exists(args.pdb_divided),args.pdb_divided
    assert os.path.exists(args.entries),args.entries
    qm = QuickModelome()
    qm.quick_structurome(xml_blast_result=args.xml_blast_result, data_dir=os.path.abspath(args.data_dir),
                         entries=args.entries, tmp_dir=args.tmp_dir,
                         pdb_divided=args.pdb_divided, max_models=args.max_models)
