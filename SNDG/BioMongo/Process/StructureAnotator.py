"""
Created on Sep 20, 2016

@author: eze
"""

import logging
import json
import os
from SNDG import mkdir

from Bio.PDB.PDBParser import PDBParser
from SNDG.BioMongo.Model.Structure import ResidueSet, ExperimentalStructure, ModeledStructure
from SNDG.Sequence.pfam import PFamProfile, ProfileNotFoundError
from SNDG.Structure.CompoundTypes import main_compound_types
from SNDG.BioMongo.Model import Structure

_log = logging.getLogger(__name__)

from SNDG.Structure.FPocket import fpocket_properties_map
eq2 = {v:k for k,v in fpocket_properties_map.items()}


class StructureAnotator(object):
    """
    classdocs
    """

    important_pfam_rs = "important_pfam"

    default_residue_sets = ["tyr", "free_try", "cyc", "free_cys", important_pfam_rs, "csa"]

    csa_prefix = "csa_"
    binding_postfix = "_binding"
    projected_postfix = "_projected"

    def __init__(self, work_dir="./", struct_path=lambda wd, strdoc: wd + strdoc.name + ".pdb"):
        """
        Constructor
        """
        self.struct_path = struct_path
        self.work_dir = work_dir


    @staticmethod
    def pocket_residue_set(pockets_json,structure_atoms):
        rss = []
        with open(pockets_json) as handle:
            pockets_dict = json.load(handle)
        for pocket_dict in pockets_dict:
            rs = ResidueSet(name="Pocket_" + str(pocket_dict["number"]), type="pocket",residues=[])
            for key, value in pocket_dict["properties"].items():
                rs[eq2[key]] = value
            residues = list(set([x.parent.parent.id + "_" + str(x.parent.id[1])
                                 for x in structure_atoms
                                 if str(x.serial_number) in pocket_dict["atoms"]]))
            rs.residues = residues
            if rs.residues:
                rss.append(rs)
        return rss


    def free_cys_tyr(self, strdoc):
        parser = PDBParser(PERMISSIVE=1, QUIET=1)

        strdoc.residue_sets = [x for x in strdoc.residue_sets if x.name not in ["free_cys", "free_tyr", "cys", "tyr"]]
        #         if not (strdoc.residue_set("free_cys") or strdoc.residue_set("free_tyr")):
        struct_path = self.struct_path(self.work_dir, strdoc)
        bp_pdb = list(parser.get_structure(strdoc.name, struct_path))[0]
        free = {"CYS": [], "TYR": []}
        _all = {"CYS": [], "TYR": []}
        codes = {"CYS": "SG", "TYR": "OH"}
        for x in bp_pdb.get_residues():
            if x.resname in codes:
                aa_code = x.parent.id + "_" + str(x.id[1])
                _all[x.resname].append(aa_code)
                neighbor_atoms = set(list(bp_pdb.get_atoms())) - set(list(x))
                if (codes[x.resname] in x) and (
                        not any(map(lambda atom: (x[codes[x.resname]] - atom) <= 3, neighbor_atoms))):
                    free[x.resname].append(aa_code)
        if free["CYS"]:
            rs = ResidueSet(name="free_cys", residues=free["CYS"])
            strdoc.residue_sets.append(rs)
        if free["TYR"]:
            rs = ResidueSet(name="free_tyr", residues=free["TYR"])
            strdoc.residue_sets.append(rs)
        if _all["CYS"]:
            rs = ResidueSet(name="cys", residues=_all["CYS"])
            strdoc.residue_sets.append(rs)
        if _all["TYR"]:
            rs = ResidueSet(name="tyr", residues=_all["TYR"])
            strdoc.residue_sets.append(rs)

    def important_pfam(self, domains, model):
        template = model.templates[0]

        important_rs = ResidueSet(name=StructureAnotator.important_pfam_rs)

        for domain in domains:

            domain_id = "_".join([domain.identifier, str(domain.location.start), str(domain.location.end)])
            domain_rs = ResidueSet(name=domain_id, type="domain",
                                   residues=[" _" + str(template.aln_pos_from_query_pos(i).q_resid)
                                             for i in range(domain.location.start, domain.location.end + 1)
                                             if i >= template.aln_query.start
                                             and i <= template.aln_query.end
                                             and i < len(template.aln_query.original_seq())
                                             ])

            model.residue_sets = [x for x in model.residue_sets if x.name != domain_id]

            if (1.0 * len(domain_rs) / len(domain)) > 0.8:
                """
                Para asegurarme que el 80% del dominio esta modelado en la estructura
                """
                model.residue_sets.append(domain_rs)

                model.residue_sets = [x for x in model.residue_sets if x.name != StructureAnotator.important_pfam_rs]
                process_important = True
                try:
                    profile = PFamProfile.create_profile(domain.identifier)
                except ProfileNotFoundError as ex:
                    _log.warn(ex)
                    process_important = False
                if process_important:
                    important = profile.important_positions()
                    if important:
                        profile_seq_map, _, aln_profile_map = profile.map_alignment_simple(domain.aln.aln_hit.seq,
                                                                                           domain.aln.aln_hit.start)
                        for i in range(len(domain.aln.aln_hit.seq)):

                            if (i in aln_profile_map) and (not domain.aln.is_gap(i)):
                                profile_pos = aln_profile_map[i]
                                if (profile_pos in important):
                                    # assert profile_pos in important
                                    seq_pos = profile_seq_map[profile_pos]
                                    important_rs.residues.append(" _" + str(seq_pos))
                        if important_rs:
                            model.residue_sets.append(important_rs)

    def calculated_annotations(self, domains, model):
        self.free_cys_tyr(model)
        try:
            self.important_pfam(domains, model)
        except IndexError as ex:
            _log.warn(ex)

    def cluster_aligned_annotation(self, model, template_aln, segment, cluster_aligned_segment, cluster_res_start):
        segment.clust_start
        segment_start = int(segment.start)
        segment_end = int(segment.end)
        segment_str = "_".join(
            [segment.pdb, segment.chain, str(segment_start).split(".")[0], str(segment.end).split(".")[0]])
        if segment_str != template_aln.aln_hit.name:
            template_eq = ExperimentalStructure.objects(name=segment.pdb).get()
            #   prot             IVAGRVSQKMAPVLRQIYDQMAEPKWVLAMGVCAS
            #   template       **IVAGAAS--MAPVLQQIYDQM--PKWVLAMGVC--**
            #   template_eq      -----AS--MAPVVQQILDQ---------------
            #   template_eq2   **IVAAAS--MAPVVQQILDQ---------------
            #   template_eq3   **IVAAAS--MAPVVQQILDQQM--PKWVLAMGVC--**

            template_res_to_cluster_res_map = lambda template_res: (
                                                                           template_res - cluster_res_start) - segment.clust_start + segment_start

            alignment_start_after_cluster_start = cluster_aligned_segment.residue_numbers()[0] >= segment.clust_start
            eq_start_residue = segment_start if alignment_start_after_cluster_start else segment_start + template_aln.aln_hit.start

            alignment_end_before_cluster_end = cluster_aligned_segment.residue_numbers()[-1] <= segment.clust_end
            eq_end_residue = segment_end if alignment_end_before_cluster_end else template_res_to_cluster_res_map(
                template_aln.hit_res_end)

            eq_aligned_segment = ResidueSet(name="aln", residues=[segment.chain + "_" + str(i + eq_start_residue)
                                                                  for i in range(len(template_aln.aln_query.seq)) if
                                                                  not template_aln.is_gap(i)])

            def map_eq_resid_to_query_resid(eq_res_id, template_aln, segment, segment_start):
                pos_in_eq = eq_res_id - segment_start
                pos_in_template = pos_in_eq + segment.clust_start
                pos_in_query = template_aln.map_pos_hit_query(pos_in_template)  # - template_aln.aln_hit.start
                return pos_in_query

            eq_aligned_csas = (template_eq.residue_set("csa").in_range(eq_start_residue,
                                                                       eq_end_residue) & eq_aligned_segment).residue_numbers()
            if eq_aligned_csas:
                aligned_csas_projected = [map_eq_resid_to_query_resid(eq_res_id, template_aln, segment, segment_start)
                                          for eq_res_id in eq_aligned_csas]
                residue_set_csa = ResidueSet(name="csa_" + segment.pdb, type="catalitic_projected",
                                             residues=["_" + str(resid) for resid in aligned_csas_projected])
                model.residue_sets.append(residue_set_csa)

            for comp_type in main_compound_types:
                binding_name = comp_type.lower() + "_binding"
                eq_aligned_binding = (template_eq.residue_set(binding_name).in_range(eq_start_residue, eq_end_residue)
                                      & eq_aligned_segment).residue_numbers()
                if eq_aligned_binding:
                    aligned_binding_projected = [
                        map_eq_resid_to_query_resid(eq_res_id, template_aln, segment, segment_start) for eq_res_id in
                        eq_aligned_binding if eq_res_id]
                    residue_set_binding = ResidueSet(name=binding_name + "_" + segment.pdb,
                                                     type=binding_name + "_projected",
                                                     residues=["_" + str(resid) for resid in aligned_binding_projected])
                    model.residue_sets.append(residue_set_binding)

    def aligned_annoations(self, model):
        model.residue_sets = [x for x in model.residue_sets
                              if not any(map(lambda rsname: x.name.startswith(rsname),
                                             ["csa"] + [y.lower() + "_binding" for y in
                                                        main_compound_types]))]

        for template_aln in model.templates:
            pdb, chain, segment_start, segment_end = template_aln.aln_hit.name.split("_")

            try:
                template = ExperimentalStructure.objects(name=pdb).get()
            except Structure.DoesNotExist as ex:
                print([ex,pdb])
                continue


            is_aligned = not ((int(segment_start) == -1) & (int(segment_end) == -1))

            if is_aligned:
                start_residue = template_aln[0].h_resid
                end_residue = template_aln[-1].h_resid
                aligned_segment = ResidueSet(name="aln", residues=[chain + "_" + str(template_aln[i].h_resid)
                                                                   for i in range(len(template_aln.aln_query.seq))
                                                                   if not template_aln.is_gap(i)])

                aligned_csas = template.residue_set("csa").in_range(start_residue, end_residue) & aligned_segment
                if len(aligned_csas):
                    aln_query = [" _" + str(template_aln.aln_pos_from_h_resid(int(x.split("_")[1])).q_resid) for x in
                                 aligned_csas.residues]
                    model.residue_sets.append(
                        ResidueSet(name="csa_" + pdb, residues=aln_query, type="catalitic_projected"))

                for comp_type in main_compound_types:
                    binding_name = comp_type.lower() + "_binding"
                    aligned_binding = template.residue_set(binding_name).in_range(start_residue,
                                                                                  end_residue) & aligned_segment

                    def query_residue(pdb_resid):
                        return " _" + str(template_aln.aln_pos_from_h_resid(pdb_resid).q_resid)

                    aligned_binding_residues = [query_residue(pdb_resid) for pdb_resid in
                                                aligned_binding.residue_numbers()]
                    if aligned_binding_residues:
                        rsbinding = ResidueSet(name=binding_name + "_" + pdb, residues=aligned_binding_residues,
                                               type=binding_name + "_projected")
                        model.residue_sets.append(rsbinding)

    #             for segment in template.clusters[0].parts:
    #                 if pdb != segment.pdb:
    #                     cluster_res_start = int(template.clusters[0].name.split("_")[2])
    #                     cluster_aligned_annotation(model, template_aln,segment,aligned_segment,cluster_res_start)

    def annotate_model(self, model, domains):
        self.calculated_annotations(domains, model)
        self.aligned_annoations(model)

    def annotate_pdb(self, model, domains):
        self.calculated_annotations(domains, model)
        self.aligned_annoations(model)



    def total(self, db, organism, query={}):
        query["organism"] = organism
        return db.structures.count(query)



    def iterator(self, db, organism, query={}):
        """
        db = MongoClient().database
        """
        query["organism"] = organism

        ids = [x["_id"] for x in db.structures.find(query, {"_id": 1})]
        for i, _id in enumerate(ids, 1):
            yield ModeledStructure.objects(id=_id).no_cache().get()
