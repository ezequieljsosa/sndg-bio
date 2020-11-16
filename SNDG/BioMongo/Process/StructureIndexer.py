# -*- coding: utf-8 -*-
"""
Created on Jun 22, 2016

@author: eze
"""
import logging
import os

import pandas as pd
from tqdm import tqdm

from SNDG.BioMongo.Model.SeqCollection import SeqColDruggabilityParam
from SNDG.BioMongo.Model.DruggabilitySearch import ProteinDruggabilitySearch, StructureDruggabilitySearch, \
    PocketDruggabilitySearch
from SNDG.BioMongo.Model.Protein import Protein
from SNDG.BioMongo.Model.SeqColDruggabilityParam import SeqColDruggabilityParamTypes
from SNDG.BioMongo.Process.StructureAnotator import StructureAnotator
from SNDG.Sequence.so import SO_TERMS
from SNDG.Structure.CompoundTypes import compound_type, main_compound_types

_log = logging.getLogger(__name__)


class StructuromeIndexer(object):
    search_params = [
        ("has_structure", "Protein gas a 3d structure", "structure",
         SeqColDruggabilityParamTypes.value, ["true", "false"], "true", "equal", "avg"),
        ("structure_type", "experimental or model", "structure",
         SeqColDruggabilityParamTypes.value, ["experimental", "model"], "experimental", "equal", "avg"),
        ("druggability",
         "Druggability score from the most druggable pocket. Druggable: druggability > 0.5 / Highly Druggable  druggability > 0.7. (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4014675/) ",
         "structure", SeqColDruggabilityParamTypes.number, None, 0.5, ">", "avg"),
        ("hydrophobicity", "Hydrophobicity of the most druggable pocket", "structure",
         SeqColDruggabilityParamTypes.number, None, None, ">", "avg"),
        ("volume", "Volume in cubic Å of the most druggable pocket",
         "structure", SeqColDruggabilityParamTypes.number, None, None, ">", "avg"),
        # ("druggability_cat","Druggability category, based on the druggability score. Druggable  (0.5 < druggability < 0.7) / Highly Druggable  (0.7 < druggability) /  ","structure",SeqColDruggabilityParamTypes.value,["poorly_druggable","druggable","highly_druggable"]),

        ("free_tyr",
         "If any of the proteins structures has a tyr with his OH oxigen atom with no surronding atoms (more than cubic Å )",
         "structure", SeqColDruggabilityParamTypes.value, ["true", "false"], None, None, None),
        ("tyr", "If any of the proteins structures has a tyr",
         "structure", SeqColDruggabilityParamTypes.value, ["true", "false"], None, None, None),
        ("free_cys",
         "If any of the proteins structures has a cys with his SH sulfur atom with no surronding atoms (more than 3 Å )",
         "structure", SeqColDruggabilityParamTypes.value, ["true", "false"], None, None, None),
        ("cys",
         "If any of the proteins structures has a cys with his SH sulfur atom with no surronding atoms (more than 3 Å )",
         "structure", SeqColDruggabilityParamTypes.value, ["true", "false"], None, None, None),
        ("csa",
         "If any of the proteins cristals or model templates, has at least one residue reported in the Catalitic Site Atlas database",
         "structure", SeqColDruggabilityParamTypes.value, ["true", "false"], None, None, None),
        ("domain_extended",
         "If any of the proteins cristals or model templates, has a domain with residues near a drug or cofactor",
         "structure", SeqColDruggabilityParamTypes.value, ["true", "false"], None, None, None),

        ("pocket_with_free_tyr",
         "If any of the proteins structures has a tyr inside a pocket with his OH oxigen atom with no surronding atoms (more than 3 Å )",
         "pocket", SeqColDruggabilityParamTypes.value, ["true", "false"], None, None, None),
        ("pocket_with_tyr", "If any of the proteins structures has a tyr inside a pocket",
         "pocket", SeqColDruggabilityParamTypes.value, ["true", "false"], None, None, None),
        ("pocket_with_free_cys",
         "If any of the proteins structures has a cys inside a pocket with his SH sulfur atom with no surronding atoms (more than 3 Å )",
         "pocket", SeqColDruggabilityParamTypes.value, ["true", "false"], None, None, None),
        ("pocket_with_cys",
         "If any of the proteins structures has a cys inside a pocket with his SH sulfur atom with no surronding atoms (more than 3 Å )",
         "pocket", SeqColDruggabilityParamTypes.value, ["true", "false"], None, None, None),
        ("pocket_with_csa",
         "If any of the proteins cristals or model templates, has at least one residue inside a pocket reported in the Catalitic Site Atlas database",
         "pocket", SeqColDruggabilityParamTypes.value, ["true", "false"], None, None, None),
        ("pocket_with_domain_extended",
         "If any of the proteins structures, has a pfam domain with residues inside a pocket near a drug or cofactor",
         "pocket", SeqColDruggabilityParamTypes.value, ["true", "false"], None, None, None),
        ("pocket_with_pfam_imp_residue",
         "If any of the proteins structures, has at least one residue, inside a pocket and a pfam domain, in contact (less than 3 Å ) with a drug or cofactor",
         "pocket", SeqColDruggabilityParamTypes.value, ["true", "false"], None, None, None)
        # TODO Estructura secundaria
    ]
    for comp_type in main_compound_types:
        comp_type_lower = comp_type.lower()
        binding_name = comp_type_lower + "_binding"
        prop = (binding_name,
                "If any of the proteins cristals or model templates, has at least one residue in contact (less than 3 Å ) with a " + comp_type_lower,
                "structure", SeqColDruggabilityParamTypes.value, ["true", "false"], "avg", "equal", "true")
        search_params.append(prop)
        prop = ("pocket_with_" + binding_name,
                "If any of the proteins cristals or model templates, has at least one pocket residue in contact (less than 3 Å ) with a" + comp_type_lower,
                "pocket", SeqColDruggabilityParamTypes.value, ["true", "false"], "avg", "equal", "true")
        search_params.append(prop)

    def __init__(self, seqcollection):
        """
        Constructor
        """
        self.collection = seqcollection
        self.user = "demo"

        if os.path.exists("/data/databases/pdb/processed/distances.tbl"):
            columns = ["pdb", "chain", "domain", "res_id", "res", "res_atom_id", "ligand_id", "ligand_name",
                       "ligand_atom", "distance"]
            df = pd.read_table("/data/databases/pdb/processed/distances.tbl", sep="\t", index_col=False, names=columns)
            self.extended_domains = {x.split("_")[0]: 1 for x in df[(df.distance < 3) & (df.domain != "NoDn") & (
                [compound_type[x] in ["DRUG", "COFACTOR"] for x in df.ligand_name])].domain}
        else:
            self.extended_domains = {}

    def is_dn_ex(self, domain):
        return domain.split(".")[0] in self.extended_domains

    def update_collection_params(self):
        for name, description, target, _type, options, defaultValue, defaultOperation, defaultGroupOperation in StructuromeIndexer.search_params:
            if not self.collection.has_druggability_param(name):
                dp = SeqColDruggabilityParam(name=name, description=description, target=target,
                                             type=_type, uploader=self.user)
                if options:
                    dp.options = options
                if defaultValue:
                    dp.defaultValue = defaultValue
                if defaultOperation:
                    dp.defaultOperation = defaultOperation
                if defaultGroupOperation:
                    dp.defaultGroupOperation = defaultGroupOperation
                self.collection.druggabilityParams.append(dp)
        self.collection.save()

    def get_protein_structures(self, protein):
        protein.search.has_structure = False
        models = protein.models()
        if len(models) > 0:
            protein.search.structure_type = "model"
        cristals = []
        for exp_struct in protein.cristals():
            cristals.append(exp_struct)
            protein.search.structure_type = "experimental"
        return cristals, models

    def annotate_aln_pocket(self, cristal, pocket, aln_pocket, ds_struct):
        ds_pocket = self.create_ds_pocket(pocket.name)
        ds_struct.pockets.append(ds_pocket)

        # Calculated props
        # Fixme, hay que actualizar los type de los residue set dominios de pdb
        dn_rss = [x for x in cristal.residue_sets if x.type == "domain" or (x.name.startswith("PF") and "_" in x.name)
                  if self.is_dn_ex(x.name)]
        ds_pocket.domain_extended = any(map(lambda dn_rs: bool(dn_rs & aln_pocket), dn_rss))
        ds_pocket.important_pfam = bool(cristal.residue_set(StructureAnotator.important_pfam_rs) & aln_pocket)
        ds_pocket.cys = bool(cristal.residue_set("cys") & aln_pocket)
        ds_pocket.tyr = bool(cristal.residue_set("tyr") & aln_pocket)
        ds_pocket.free_cys = bool(cristal.residue_set("free_cys") & aln_pocket)
        ds_pocket.free_tyr = bool(cristal.residue_set("free_tyr") & aln_pocket)

        ds_pocket.druggability = pocket.get_druggability_score()
        ds_pocket.hydrophobicity = pocket.get_hydrophobicity_score()
        ds_pocket.volume = pocket.get_volume()

        # Aligned Props
        for comp_type in main_compound_types:
            comp_type_lower = comp_type.lower()
            binding_name = comp_type_lower + "_binding"
            ds_pocket[comp_type_lower] = bool(cristal.residue_set(binding_name) & aln_pocket)
        ds_pocket.csa = bool(cristal.residue_set("csa") & aln_pocket) and (ds_pocket.druggability > 0.5)
        return ds_pocket

    def annotate_modeled_pocket(self, model, pocket, ds_struct):
        ds_pocket = self.create_ds_pocket(pocket.name)

        ds_struct.pockets.append(ds_pocket)

        dn_rss = [x for x in model.residue_sets if x.type == "domain" and self.is_dn_ex(x.name)]
        ds_pocket.domain_extended = any(map(lambda dn_rs: bool(pocket & dn_rs), dn_rss))

        ds_pocket.important_pfam = bool(model.residue_set(StructureAnotator.important_pfam_rs) & pocket)
        ds_pocket.cys = bool(model.residue_set("cys") & pocket)
        ds_pocket.tyr = bool(model.residue_set("tyr") & pocket)
        ds_pocket.free_cys = bool(model.residue_set("free_cys") & pocket)
        ds_pocket.free_tyr = bool(model.residue_set("free_tyr") & pocket)

        ds_pocket.druggability = pocket.get_druggability_score()
        ds_pocket.hydrophobicity = pocket.get_hydrophobicity_score()
        ds_pocket.volume = pocket.get_volume()

        for template in model.template_structures():
            self.annotate_pocket_with_template(pocket, model, template.name, ds_pocket)

        return ds_pocket

    def annotate_with_pocket(self, ds_prot, ds_struct, ds_pocket):

        ds_prot.pocket_with_free_tyr = ds_prot.pocket_with_free_tyr | ds_pocket.free_tyr
        ds_prot.pocket_with_free_cys = ds_prot.pocket_with_free_cys | ds_pocket.free_cys
        ds_prot.pocket_with_tyr = ds_prot.pocket_with_tyr | ds_pocket.tyr
        ds_prot.pocket_with_cys = ds_prot.pocket_with_cys | ds_pocket.cys
        ds_prot.pocket_with_pfam_imp_residue = ds_prot.pocket_with_pfam_imp_residue | ds_pocket.important_pfam
        ds_prot.pocket_with_csa = ds_prot.pocket_with_csa | ds_pocket.csa
        ds_prot.pocket_with_domain_extended = ds_prot.pocket_with_domain_extended | ds_pocket.domain_extended

        ds_struct.pocket_with_free_tyr = ds_struct.pocket_with_free_tyr | ds_pocket.free_tyr
        ds_struct.pocket_with_free_cys = ds_struct.pocket_with_free_cys | ds_pocket.free_cys
        ds_struct.pocket_with_tyr = ds_struct.pocket_with_tyr | ds_pocket.tyr
        ds_struct.pocket_with_cys = ds_struct.pocket_with_cys | ds_pocket.cys
        ds_struct.pocket_with_pfam_imp_residue = ds_struct.pocket_with_pfam_imp_residue | ds_pocket.important_pfam
        ds_struct.pocket_with_csa = ds_struct.pocket_with_csa | ds_pocket.csa
        ds_struct.pocket_with_domain_extended = ds_struct.pocket_with_domain_extended | ds_pocket.domain_extended

        for comp_type in main_compound_types:
            comp_type_lower = comp_type.lower()
            binding_name = comp_type_lower + "_binding"
            ds_prot["pocket_with_" + binding_name] = ds_prot["pocket_with_" + binding_name] | ds_pocket[binding_name]
            ds_struct["pocket_with_" + binding_name] = ds_struct["pocket_with_" + binding_name] | ds_pocket[
                binding_name]

        if ds_prot.druggability < ds_pocket.druggability:
            ds_prot.druggability = ds_pocket.druggability
            ds_prot.hydrophobicity = ds_pocket.hydrophobicity
            ds_prot.volume = ds_pocket.volume

        if ds_struct.druggability < ds_pocket.druggability:
            ds_struct.druggability = ds_pocket.druggability
            ds_struct.hydrophobicity = ds_pocket.hydrophobicity
            ds_struct.volume = ds_pocket.volume

    def process_cristal(self, protein, cristal):
        i = 0
        for f in protein.features:

            if f.type == SO_TERMS["polypeptide_structural_motif"] and cristal.name in f.identifier:
                i = i + 1
                if i > 30:
                    break
                self.process_cristal_feature(protein, cristal, f)

    def process_cristal_feature(self, protein, cristal, feature):
        hit_arr = feature.aln.aln_hit.name.split("_")
        if (len(hit_arr) > 4) and hit_arr[4].startswith("PF"):
            _, chain, _, _, _, dnstart, _ = feature.aln.aln_hit.name.split("_")
            offset = int(dnstart)
        else:
            _, chain = feature.aln.aln_hit.name.split("_")[0:2]
            offset = 0

        aln_residue_set = feature.aln.residue_set_aln(cristal, chain, offset=offset)
        if aln_residue_set:
            ds_prot = protein.search
            ds_struct = StructureDruggabilitySearch(structure=feature.aln.aln_hit.name, druggability=0)
            self.initDrugabilitySearch(ds_struct)
            ds_prot.structures.append(ds_struct)

            aligned_pockets = [(p, p & aln_residue_set) for p in cristal.pockets if p & aln_residue_set]

            self.annotate_cristal(cristal, ds_struct, aln_residue_set, ds_prot)
            self.annotate_protein_with_cristal(protein, ds_struct)

            for pocket, aln_pocket in aligned_pockets:
                ds_pocket = self.annotate_aln_pocket(cristal, pocket, aln_pocket, ds_struct)
                self.annotate_with_pocket(ds_prot, ds_struct, ds_pocket)
        else:
            protein.features = [f for f in protein.features if f != feature]

    def process_model(self, protein, model):

        ds_prot = protein.search
        ds_struct = StructureDruggabilitySearch(structure=model.name, druggability=0)
        self.initDrugabilitySearch(ds_struct)
        ds_prot.structures.append(ds_struct)

        if hasattr(model, "templates"):
            protein.keywords.append(model.templates[0].aln_hit.name.lower())

        self.annotate_structure(ds_prot, model, ds_struct)

        for template in model.templates:
            template_name = template.aln_hit.name.split("_")[0]
            self.annotate_with_template(protein, model, template_name, ds_struct)

        for pocket in model.pockets:
            ds_pocket = self.annotate_modeled_pocket(model, pocket, ds_struct)
            self.annotate_with_pocket(ds_prot, ds_struct, ds_pocket)

    def annotate_pocket_with_template(self, pocket, model, template_name, ds_pocket):

        for comp_type in main_compound_types:
            comp_type_lower = comp_type.lower()
            search_comp = comp_type_lower + "_binding"
            binding_name = search_comp + "_" + template_name

            pocket_has_comp = bool(model.residue_set(binding_name) & pocket)
            if hasattr(ds_pocket, search_comp):
                ds_pocket[search_comp] = ds_pocket[search_comp] | pocket_has_comp
            else:
                ds_pocket[search_comp] = pocket_has_comp
        has_csa = (bool(model.residue_set("csa_" + template_name) & pocket) and (pocket.druggability_score > 0.5))
        if hasattr(ds_pocket, "csa"):
            ds_pocket.csa = ds_pocket.csa | has_csa
        else:
            ds_pocket.csa = has_csa

    def annotate_cristal(self, cristal, ds_struct, aln_residue_set, ds_prot):

        for comp_type in main_compound_types:
            comp_type_lower = comp_type.lower()
            binding_name = comp_type_lower + "_binding"
            ds_struct[binding_name] = bool(cristal.residue_set(binding_name) & aln_residue_set)

        ds_struct.csa = bool(cristal.residue_set("csa") & aln_residue_set)
        self.annotate_structure(ds_prot, cristal, ds_struct)

    def annotate_structure(self, ds_prot, structure, ds_struct):

        ds_struct.domain_extended = any(
            map(lambda rs: self.is_dn_ex(rs.name), structure.residue_sets_for_type("domain")))

        ds_struct.important_pfam = bool(structure.residue_set(StructureAnotator.important_pfam_rs))
        ds_struct.cys = bool(structure.residue_set("cys"))
        ds_struct.tyr = bool(structure.residue_set("tyr"))
        ds_struct.free_cys = bool(structure.residue_set("free_cys"))
        ds_struct.free_tyr = bool(structure.residue_set("free_tyr"))

        ds_prot.domain_extended = ds_prot.domain_extended | ds_struct.domain_extended
        ds_prot.free_tyr = ds_prot.free_tyr | ds_struct.free_tyr
        ds_prot.free_cys = ds_prot.free_cys | ds_struct.free_cys
        ds_prot.tyr = ds_prot.tyr | ds_struct.tyr
        ds_prot.cys = ds_prot.cys | ds_struct.cys

        for ds_pocket in ds_struct.pockets:
            if ds_struct.druggability < ds_pocket.druggability:
                ds_struct.druggability = ds_pocket.druggability
                ds_struct.hydrophobicity = ds_pocket.hydrophobicity
                ds_struct.volume = ds_pocket.volume

    def annotate_protein_with_cristal(self, protein, ds_struct):
        ds_prot = protein.search
        ds_prot.csa = ds_prot.csa | ds_struct.csa

        ds_prot.free_tyr = ds_prot.free_tyr | ds_struct.free_tyr
        ds_prot.free_cys = ds_prot.free_cys | ds_struct.free_cys
        ds_prot.tyr = ds_prot.tyr | ds_struct.tyr
        ds_prot.cys = ds_prot.cys | ds_struct.cys

        prot_extended_domains = [y.identifier for y in protein.domains() if
                                 y.identifier.split(".")[0] in self.extended_domains]
        ds_prot.domain_extended = bool(prot_extended_domains)
        for comp_type in main_compound_types:
            comp_type_lower = comp_type.lower()
            binding_name = comp_type_lower + "_binding"
            ds_prot[binding_name] = ds_prot[binding_name] | ds_struct[binding_name]

    def annotate_with_template(self, protein, structure, template_name, ds_struct):
        ds_prot = protein.search

        for comp_type in main_compound_types:
            comp_type_lower = comp_type.lower()
            binding_name = comp_type_lower + "_binding"

            if hasattr(ds_struct, binding_name):
                ds_struct[binding_name] = ds_struct[binding_name] | bool(
                    structure.residue_set(binding_name + "_" + template_name))
            else:
                ds_struct[binding_name] = bool(structure.residue_set(binding_name + "_" + template_name))

            ds_prot[binding_name] = ds_prot[binding_name] | ds_struct[binding_name]
        ds_prot.csa = ds_prot.csa | bool(structure.residue_set("csa_" + template_name))
        ds_prot.free_tyr = ds_prot.free_tyr | ds_struct.free_tyr
        ds_prot.free_cys = ds_prot.free_cys | ds_struct.free_cys
        ds_prot.tyr = ds_prot.tyr | ds_struct.tyr
        ds_prot.cys = ds_prot.cys | ds_struct.cys

    def annotate_protein_with_structures(self, cristals, models, protein):

        for cristal in cristals:
            try:
                self.process_cristal(protein, cristal)
            except Exception as ex:
                _log.warn(ex)
        for model in models:
            try:
                self.process_model(protein, model)
            except Exception as ex:
                _log.warn(ex)

        protein.search.has_structure = len(protein.structures()) > 0
        if protein.search.has_structure:
            protein.keywords.append("has_structure")
            protein.keywords = [x for x in protein.keywords if
                                x not in ["poorly_druggable", "druggable", "highly_druggable"]]
        else:
            protein.search.druggability = 0

        drugability = "non_druggable"
        if (protein.search.druggability >= 0.2) and (protein.search.druggability < 0.5):
            drugability = "poorly_druggable"
        elif (protein.search.druggability >= 0.5) and (protein.search.druggability < 0.7):
            drugability = "druggable"
        elif protein.search.druggability > 0.7:
            drugability = "highly_druggable"

        protein.search.druggability_cat = drugability
        protein.keywords = list(set(protein.keywords + [drugability]))

    def create_ds_pocket(self, pocket_name):
        ds = PocketDruggabilitySearch(pocket=pocket_name)
        ds.druggability = 0

        for comp_type in main_compound_types:
            comp_type_lower = comp_type.lower()
            binding_name = comp_type_lower + "_binding"
            ds[binding_name] = False
            # ocket_with_drug_binding
            ds["pocket_with_" + binding_name] = False

        ds.domain_extended = False
        ds.csa = False

        ds.free_tyr = False
        ds.free_cys = False
        ds.tyr = False
        ds.cys = False

        return ds

    def initDrugabilitySearch(self, ds):

        ds.pocket_with_free_tyr = False
        ds.pocket_with_free_cys = False
        ds.pocket_with_tyr = False
        ds.pocket_with_cys = False
        ds.pocket_with_pfam_imp_residue = False
        ds.pocket_with_csa = False
        ds.pocket_with_domain_extended = False
        ds.druggability = 0

        for comp_type in main_compound_types:
            comp_type_lower = comp_type.lower()
            binding_name = comp_type_lower + "_binding"
            ds[binding_name] = False
            # ocket_with_drug_binding
            ds["pocket_with_" + binding_name] = False

        ds.domain_extended = False
        ds.csa = False

        ds.free_tyr = False
        ds.free_cys = False
        ds.tyr = False
        ds.cys = False

        return ds

    def build_index(self, query={}, error_output="/tmp/index_struct.log"):
        self.update_collection_params()
        if not query:
            proteins = Protein.objects(organism=self.collection.name).no_cache().timeout(False)
            prot_count = Protein.objects(organism=self.collection.name).count()
        else:
            query["organism"] = self.collection.name
            proteins = Protein.objects(__raw__=query).no_cache().timeout(False)
            prot_count = Protein.objects(__raw__=query).count()
        with tqdm(proteins, total=prot_count) as pbar:
            for protein in pbar:

                pbar.set_description(protein.name)
                try:
                    if not protein.search:
                        protein.search = ProteinDruggabilitySearch()

                    self.initDrugabilitySearch(protein.search)

                    cristals, models = self.get_protein_structures(protein)

                    if len(cristals + models) > 0:
                        self.annotate_protein_with_structures(cristals, models, protein)

                    else:
                        protein.keywords = [x for x in protein.keywords if
                                            x not in ["poorly_druggable", "druggable", "highly_druggable",
                                                      "has_structure"]]
                        protein.keywords.append("non_druggable")
                        for x in StructuromeIndexer.search_params:
                            protein.search[x[0]] = None
                    protein.save()
                except Exception as ex:
                    error_line = protein.name + "," + protein.organism + "," + str(protein.id) + "," + str(ex)
                    _log.warn(error_line)
                    with open(error_output, "a") as h:
                        h.write(error_line)
