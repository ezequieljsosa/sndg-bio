"""

"""
import glob
import os
from tempfile import mkdtemp

import Bio.SearchIO as bpsio
from SNDG import execute
from SNDG.Structure.Modelome import Modelome

from modeller.automodel import assess, refine  # @UnresolvedImport
REFINEMENT = refine.very_slow
ASSESMENTS = [assess.GA341, assess.DOPE]  # @UndefinedVariable
MODELS_TO_GENERATE = 1
import logging
_log = logging.getLogger(__name__)

class PsiProfile:
    @staticmethod
    def build_profile(seq_fasta, database, iterations, pssm_file, cpu):
        execute(
            "psiblast -query {input} -db {database} -num_threads {cpu} -out_pssm {output} -evalue 0.0001 -num_iterations {iterations} > {cmd_out}",
            output=pssm_file, iterations=iterations,
            database=database, input=seq_fasta, cpu=cpu, cmd_out=pssm_file + ".out")

    @staticmethod
    def profile_search(seq_id, database, pssm_file, search_result, cpu):
        execute(
            "psiblast -db {database} -in_pssm {input} -num_threads {cpu}  -evalue 0.001  -outfmt 5 -out {output} > {cmd_out} ",
            output=search_result,
            database=database, input=pssm_file, cpu=cpu, cmd_out=search_result + ".out")

    @staticmethod
    def create_psi_model(seq_id, seq_fasta, profile_database, pssm_build_iterations, pdb_seqres, work_dir, cpus,
                         refinement=REFINEMENT, models_to_generate=MODELS_TO_GENERATE, assessments=ASSESMENTS,
                         entries={},
                         pdb_divided="/data/databases/pdb/divided/", tmp_dir=None,skip_profile=False,skip_quality=False):
        if not tmp_dir:
            tmp_dir = mkdtemp("structurome_")
        search_result = work_dir + "/profile_search.xml"

        pssm_file = work_dir + "/profile.pssm"

        if skip_profile and not os.path.exists(pssm_file):
            _log.warning("no pssm for " + seq_id)
            return

        if not os.path.exists(pssm_file):
            PsiProfile.build_profile(seq_fasta, profile_database, pssm_build_iterations, pssm_file, cpus)

        PsiProfile.profile_search(seq_id, pdb_seqres, pssm_file, search_result, cpus)
        search_result = bpsio.parse(search_result, "blast-xml")
        hsps = []
        for query in search_result:
            for hit in list(query):
                for hsp in hit:
                    if hsp.evalue < 10 ** -5:
                        hsps.append(hsp)
        result = Modelome.model_hsps(seq_id, work_dir, hsps, refinement, models_to_generate, assessments, entries, pdb_divided,
                            tmp_dir,max_models=3)

        if not skip_quality:
            for models in result["models"].values():
                for model_path in models:
                    if not os.path.exists(model_path + ".json"):
                        assessment = QMean.assesment(model_path)
                        with open(model_path + ".json", "w") as h:
                            json.dump(assessment, h)

    def protein_models(self,work_dir,protein_name):
        return glob.glob(work_dir + "/" + protein_name + "/" + protein_name + "/*")