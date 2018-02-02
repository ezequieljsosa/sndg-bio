"""

"""

import os
import json
from tempfile import mkdtemp
import Bio.SearchIO as bpsio

from SNDG import execute, Struct
from SNDG.Structure.Modeller import Modeller
from SNDG.Structure.QMean import QMean
from SNDG.Structure.ChainSplitter import ChainSplitter


def build_profile(seq_fasta, database, iterations, pssm_file, cpu):
    execute(
        "psiblast -query {input} -db {database} -num_threads {cpu} -out_pssm {output} -evalue 0.0001 -num_iterations {iterations} > {cmd_out}",
        output=pssm_file, iterations=iterations,
        database=database, input=seq_fasta, cpu=cpu, cmd_out=pssm_file + ".out")


def profile_search(seq_id, database, pssm_file, search_result, cpu):
    execute(
        "psiblast -db {database} -in_pssm {input} -num_threads {cpu}  -evalue 0.001  -outfmt 5 -out {output} > {cmd_out} ",
        output=search_result,
        database=database, input=pssm_file, cpu=cpu, cmd_out=search_result + ".out")


from modeller.automodel import assess, refine  # @UnresolvedImport

REFINEMENT = refine.very_slow
ASSESMENTS = [assess.GA341, assess.DOPE]  # @UndefinedVariable
MODELS_TO_GENERATE = 1


def create_psi_model(seq_id, seq_fasta, profile_database, pssm_build_iterations, pdb_seqres, work_dir, cpus,
                     refinement=REFINEMENT, models_to_generate=MODELS_TO_GENERATE, assessments=ASSESMENTS, entries={},
                     pdb_divided="/data/databases/pdb/divided/", tmp_dir=None):
    if not tmp_dir:
        tmp_dir = mkdtemp("structurome_")
    search_result = work_dir + "/profile_search.xml"
    pssm_file = work_dir + "/profile.pssm"
    if not os.path.exists(pssm_file):
        build_profile(seq_fasta, profile_database, pssm_build_iterations, pssm_file, cpus)

    profile_search(seq_id, pdb_seqres, pssm_file, search_result, cpus)
    search_result = bpsio.parse(search_result, "blast-xml")
    hsps = []
    for query in search_result:
        for hit in list(query):
            for hsp in hit:
                if hsp.evalue < 10 ** -5:
                    hsps.append(hsp)
    model_hsps(seq_id,work_dir, hsps, refinement, models_to_generate, assessments, entries, pdb_divided, tmp_dir)


def model_hsps(seq_id,work_dir, hsps, refinement=REFINEMENT, models_to_generate=MODELS_TO_GENERATE, assessments=ASSESMENTS,
               entries={}, pdb_divided="/data/databases/pdb/divided/", tmp_dir=None):
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
    for aligment in alns[0:3]:
        # pdb,aligment = pdb_alignment
        pdb, chain, _, _ = aligment.aln_hit.name.split("_")
        if not os.path.exists(modeler.pdb_path(seq_id + "_" + aligment.aln_hit.name, seq_id)):
            base_model_path = pdb_divided + pdb[1:3] + "/pdb" + pdb + ".ent"
            ChainSplitter(tmp_dir).make_pdb(base_model_path, pdb, chain, overwrite=True)
            models = modeler.create_model(seq_id + "_" + aligment.aln_hit.name, aligment)
        else:
            models = [modeler.pdb_path(seq_id + "_" + aligment.aln_hit.name, seq_id, idx) for idx in
                      range(1, models_to_generate + 1)]

        for model_path in models:
            if not os.path.exists(model_path + ".json"):
                assessment = QMean.assesment(model_path)
                with open(model_path + ".json", "w") as h:
                    json.dump(assessment, h)

