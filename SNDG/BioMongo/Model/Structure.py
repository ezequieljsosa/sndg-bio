'''
Created on Oct 9, 2015

@author: eze
'''

import itertools
import logging

from Bio.SeqUtils import seq1
from mongoengine.document import Document, EmbeddedDocument, \
    DynamicEmbeddedDocument
from mongoengine.errors import DoesNotExist
from mongoengine.fields import StringField, ListField, EmbeddedDocumentField, \
    ReferenceField, FloatField, IntField,DynamicField,DictField

from SNDG.BioMongo.Model import BioProperties
from SNDG.BioMongo.Model.Alignment import SimpleAlignment
from SNDG.BioMongo.Model.SeqCollection import SeqCollection
from SNDG.BioMongo.Model.Protein import Protein
from SNDG.Sequence.so import SO_TERMS
from SNDG.Structure.CompoundTypes import compound_type
from SNDG.BioMongo.Model import Cluster
from SNDG.BioMongo.Model.ResidueAln import ResidueAln #ignore unused




_log = logging.getLogger(__name__)

PDB_LIGAND_METALS = [x for x, y in compound_type.items() if y == "METAL"]


def cristals(self, raise_error_if_missing=False):
    def get_cristal(identifier):
        try:
            return ExperimentalStructure.objects(name=identifier.split("_")[0]).no_cache().get()
        except DoesNotExist:
            if raise_error_if_missing:
                raise
            _log.warn("Structure %s does not exists" % identifier)
            return None

    return [x for x in [get_cristal(feature.identifier)
                        for feature in self.features if feature.type == SO_TERMS["polypeptide_structural_motif"]] if x]


setattr(Protein, 'cristals', cristals)


def models(self):
    for alias in self.alias:
        models = list(ModeledStructure.objects(organism=self.organism, templates__aln_query__name=alias).no_cache())
        if models:
            return models
    return []


setattr(Protein, 'models', models)


def structures(self):
    return set(self.cristals() + self.models())


setattr(Protein, 'structures', structures)


def residue_set_aln(self, structure, chain_name, offset=0):
    # TODO deberia tener menos parametros??
    rs = ResidueSet(name="aln_" + structure.name + "_" + chain_name)
    chain = structure.chain(chain_name)
    delta = 0
    mol = chain.residues[self.aln_hit.start + offset]
    if seq1(mol.compound).lower() != self.aln_hit.seq.replace("-", "")[0].lower():
        delta = -1

    error = 0
    for i, aa in enumerate(self.aln_hit.seq.replace("-", "")):
        try:
            mol = chain.residues[self.aln_hit.start + i + delta + offset]
            if seq1(mol.compound).lower() != aa.lower():
                error += 1
                if error > 10:
                    raise Exception("too many mismaches")
            #                 assert seq1( mol.compound).lower() == aa.lower(), (mol.compound, aa)
            residue = chain_name + "_" + str(mol.resid)
            rs.residues.append(residue)
        except Exception as ex:
            _log.error("pdb %s mal alineado con residuos: %s" % (self.aln_hit.name, ex))
            return ResidueSet(name="aln_" + structure.name + "_" + chain_name)

    return rs


setattr(SimpleAlignment, 'residue_set_aln', residue_set_aln)


class ResidueSet(DynamicEmbeddedDocument):
    meta = {'allow_inheritance': True}
    name = StringField(required=True)
    type = StringField(required=False)
    #     properties = EmbeddedDocumentField(BioProperties)
    residues = ListField(StringField(), default=[])

    def __init__(self, **kwargs):
        '''
        '''

        super(EmbeddedDocument, self).__init__(**kwargs)

        self._instance = None


    def __and__(self, other):
        assert other != None, "cant & with None"
        return ResidueSet(name=self.name + "&" + other.name,
                          residues=list(
                              set([x.strip() for x in self.residues]) & set([x.strip() for x in other.residues])))

    def __nonzero__(self):
        return len(self._data["residues"]) > 0

    def __repr__(self):
        return "RS(" + self.name + " type=" + str(self.type) + ", count: " + str(len(self.residues)) + ")"

    def __str__(self):
        return self.__repr__()

    def __len__(self):
        return len(self._data["residues"])

    def residue_numbers(self):
        return [int(x.split("_")[1]) for x in self.residues]

    def in_range(self, start, end):
        return ResidueSet(name=self.name, residues=[x for x in self.residues if start <= int(x.split("_")[1]) <= end])

    def get_druggability_score(self):
        if hasattr(self, "druggability_score"):
            return self.druggability_score
        else:
            return self['Druggability Score']

    def get_hydrophobicity_score(self):
        if hasattr(self, "hydrophobicity_score"):
            return self.druggability_score
        else:
            return self['Hydrophobicity score']

    def get_volume(self):
        if hasattr(self, "volume"):
            return self.volume
        else:
            return self['Volume']


class Molecule(EmbeddedDocument):
    chain = StringField()
    compound_type = StringField(max_length=15)  # , choices=COMPOUND_TYPE )
    compound = StringField(max_length=3)
    resid = IntField(required=True)
#     sasas = MapField(FloatField())


class Chain(EmbeddedDocument):
    name = StringField(required=True)
    segments = ListField(ListField(IntField()))
    aln = EmbeddedDocumentField(SimpleAlignment)
    residues = ListField(EmbeddedDocumentField(Molecule), default=[])


class StructureQuality(EmbeddedDocument):
    name = StringField(required=True)
    value = FloatField(required=True)
    properties = EmbeddedDocumentField(BioProperties)


class Structure(Document):
    meta = {'allow_inheritance': True, 'collection': "structures",
            'index_cls': False,
            'indexes': [
                "name",
                {"fields": ["organism", "name"]},
                {"fields": ["organism", "chains.aln_query.name"]},
                {"fields": ["chains.aln_query.name"]},
                {"fields": ["seq_collection_id", "name"]},
                {"fields": ["seq_collection_name", "name"]},
                {"fields": ["_cls", "name"]},
                {"fields": ["seq_collection_id", "chains.aln.aln_query"]}],
            'db_alias': 'pdb'
            }
    '''
    classdocs
    name, description, organism , seq_collection_id , 
    chains
    residue_sets    
    file_paths
    ligands
    properties    
    qualities
    '''
    name = StringField(required=True)
    description = StringField()
    organism = StringField()
    seq_collection_name = StringField()
    seq_collection_id = ReferenceField(SeqCollection)
    chains = ListField(EmbeddedDocumentField(Chain))
    residue_sets = ListField(EmbeddedDocumentField(ResidueSet))
    pockets = ListField(EmbeddedDocumentField(ResidueSet))
    ligands = ListField(EmbeddedDocumentField(Molecule), default=[])
    properties = EmbeddedDocumentField(BioProperties)
    keywords = ListField(StringField(),default=[])
    qualities = ListField(EmbeddedDocumentField(StructureQuality))
    sndg_index = DynamicField(required=False)

    def chain(self, chain_name):
        return [x for x in self.chains if x.name == chain_name][0]

    def druggability(self):
        druggabilities = [x.druggability_score for x in self.pockets]
        if druggabilities:
            return max(druggabilities)
        return 0

    def new_residue_set(self, compound_type):
        rs = ResidueSet(name=compound_type, residues=[])
        for chain in self.chains:
            for residue in chain.residues:
                if residue.compound_type == compound_type:
                    rs.residues.append(residue.chain + "_" + str(residue.resid))
        return rs

    def residue_set(self, rs_name):
        rss = [x for x in self.residue_sets if x.name == rs_name]
        if rss:
            return rss[0]
        else:
            return ResidueSet(name=rs_name)

    def residue_sets_for_type(self, rs_type):
        return [x for x in self.residue_sets if x.type == rs_type]

    def has_residue_set(self, rs_name):
        rss = [x for x in self.residue_sets if x.name == rs_name]
        return True if len(rss) else False

    def isResidueFromPocket(self, chain, res_id):
        for p in self.pockets:
            if chain + "_" + str(res_id) in p.residues:
                return True
        return False

    def get_pocket_from_residue(self, chain, res_id):
        for p in self.pockets:
            if chain + "_" + str(res_id) in p.residues:
                return p
        raise Exception("not found: %s in %s pockets" % (chain + "_" + str(res_id), self.name))

    def all_intersections(self):
        residue_sets = self.residue_sets + [self.pockets]
        sets_of_residues_map = {}
        for residue_set in residue_sets:
            for residue in residue_set.residues:
                if residue not in sets_of_residues_map:
                    sets_of_residues_map[residue] = []
                sets_of_residues_map[residue].append(residue_set.name)
        intersections = {}
        for residue, sets in sets_of_residues_map.items():
            subsets = [tuple(sorted(sets))] if len(sets) > 1 else []
            for subset in range(len(sets))[2:]:
                subsets = subsets + [tuple(sorted(x)) for x in itertools.combinations(sets, subset)]
            subsets = set(subsets)

            for subset in subsets:
                if subset not in intersections:
                    intersections[subset] = []
                intersections[subset].append(residue)
        return intersections

    def has_metal(self):
        return any([x.compound_type in PDB_LIGAND_METALS for x in self.ligands])

    def quality(self, metric_name):
        for x in self.qualities:
            if x.name == metric_name:
                return x.value
        return None


class ExperimentalStructure(Structure):
    quaternary = StringField()
    resolution = FloatField()
    experiment = StringField()
    clusters = ListField(EmbeddedDocumentField(Cluster))
    tax = IntField(required=False)


    def cluster(self, cluster_name):
        rss = [x for x in self.clusters if x.name == cluster_name]
        if rss:
            return rss[0]
        else:
            return Cluster(name=cluster_name)

    def file_path(self):
        return '/data/pdb/divided/' + self.name[1:3] + '/pdb' + self.name + '.ent'

    def __str__(self):
        return "ExperimentalStructure(" + self.name + ")"

    def __repr__(self):
        return self.__str__()


class ModeledStructure(Structure):
    meta = {'allow_inheritance': True,
            'indexes': [
                {"fields": ["organism", "templates.aln_query.name"]}
            ]
            }

    templates = ListField(EmbeddedDocumentField(SimpleAlignment), default=[])
    pipeline = StringField(required=True)

    def __str__(self):
        return "ModeledStructure(" + self.name + ", pipe=" + self.pipeline + ")"

    def file_path(self):
        return '/data/organismos/' + self.organism + '/estructura/' + self.pipeline + "/modelos/" + self.name + '.pdb'

    def template_structures(self):
        templates = []
        for x in self.templates:
            try:
                struct = ExperimentalStructure.objects(name=x.aln_hit.name.split("_")[0]).get()
                templates.append(struct)
            except DoesNotExist as ex:
                print  ([ex, x.aln_hit.name.split("_")[0]])
        return templates

