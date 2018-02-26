"""

"""

import glob
import logging
from collections import defaultdict

from tqdm import tqdm

import Bio.SearchIO as bpsio
from SNDG.Structure.Modelome import Modelome

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
        self.hsp_dict = defaultdict(lambda : {})


    def load_hsp_dict(self,xml_blast_result):

        for query in bpsio.parse(xml_blast_result, "blast-xml"):
            for hit in query:
                if list(hit):
                    hsp = list(hit)[0]
                    self.hsp_dict[query.id][hsp.hit.id] = hsp


    @staticmethod
    def protein_models( work_dir, protein_name):
        return glob.glob(work_dir + "/" + protein_name + "/" + "/*")

    def aln_hsp(self, base_path,prot,template):
        return self.hsp_dict[prot][template]

    def quick_structurome(self, xml_blast_result, data_dir, entries, tmp_dir="/tmp/chain_PDBs"):

        con_pdb = []
        good_model = defaultdict(lambda: [])

        def identity(hsp):
            return 1.0 * hsp.ident_num / hsp.aln_span

        for query in bpsio.parse(xml_blast_result, "blast-xml"):
            for hit in query:
                if list(hit):
                    hsp = list(hit)[0]
                    if identity(hsp) >= 0.9:
                        con_pdb.append(hsp.query.id)
                    elif identity >= 0.6:
                        good_model[hsp.query.id].append(hsp)

        tuplas = good_model.items()

        with tqdm(tuplas) as pbar:
            for seq, hsps in pbar:
                try:
                    from SNDG.Structure.Modelome import Modelome
                    Modelome.model_hsps(seq, data_dir, hsps, entries=entries, tmp_dir=tmp_dir)
                except Exception as ex:
                    _log.exception(ex)
