import copy

from mongoengine import Document, EmbeddedDocument, DynamicEmbeddedDocument
from mongoengine.fields import StringField, ListField, EmbeddedDocumentField, ObjectIdField, IntField, \
    DynamicField, ReferenceField, BinaryField

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from SNDG.BioMongo.Model.DruggabilitySearch import ProteinDruggabilitySearch
from SNDG.BioMongo.Model.Feature import Feature
from SNDG.BioMongo.Model.SeqCollection import SeqCollection


class NotFoundException(Exception):
    pass


class Size(EmbeddedDocument):
    unit = StringField(max_length=10)
    len = IntField(required=True)

class EmbeddedSNDGIndex(DynamicEmbeddedDocument):
        tax = ListField(StringField())


class BioProperty(DynamicEmbeddedDocument):
    '''
    '''
    #     field_list = ["_type","property","value","description","url","source"]

    _type = StringField(required=True)
    property = StringField()
    value = DynamicField()
    description = StringField()
    url = StringField()
    source = StringField()

    def __init__(self, **kwargs):
        '''
        '''
        super(DynamicEmbeddedDocument, self).__init__(**kwargs)

    def __str__(self):
        return "BioProperty(_type={_type},property={property}, value={value}) ".format(
            _type=str(self._type), property=str(len(self.property)), value=self.value)

    def __repr__(self):
        return self.__str__()


class Sequence(Document):
    '''
    Represents a collections' attributes, not the list itself.
    * collection_type: genome(SO:0001026),peptide_collection(SO:0001501),
             contig_collection(SO:0001462), variant_collection(SO:0001507)
    * status: chromosome(SO:0000340),scaffold(SO:0000148) or contig("SO:0000149). See SO_TERMS_STATUS
    * organism: the sequences from the collection belongs to this organism
    * name: Step name
    '''

    meta = {'allow_inheritance': True, 'abstract': True}

    name = StringField(required=True)
    description = StringField(required=False)
    seq = StringField(default="", required=True)
    features = ListField(EmbeddedDocumentField(Feature))
    properties = ListField(EmbeddedDocumentField(BioProperty))
    ontologies = ListField(StringField())
    size = EmbeddedDocumentField(Size)
    status = StringField(default="predicted")
    keywords = ListField(StringField())
    dbxrefs = ListField(StringField())
    auth = StringField()

    search = EmbeddedDocumentField(ProteinDruggabilitySearch)
    sndg_index = EmbeddedDocumentField(EmbeddedSNDGIndex,required=False)

    def __init__(self, **kwargs):
        '''
        '''
        super(Document, self).__init__(**kwargs)
        self._seq_init()

    def _seq_init(self):
        self._seqrecord = SeqRecord(Seq(self.seq), str(self.id), self.name)
        self._seqrecord.description = self.description
        self.size = Size(unit="?", len=len(self.seq))
        self.features = self.features if self.features else []

    def add_feature(self, feature):
        self.features.append(feature)

    def _feature_name(self, feature):
        return "_".join([self.name, feature.identifier, str(feature.location.start), str(feature.location.end)])

    def _feature_id(self, feature):
        if feature._id:
            return feature._id
        else:
            return self._feature_name(feature)

    def seqs_from_features(self, **kwargs):
        search_fn = lambda f, x: True
        if kwargs:
            search_fn = lambda feature, search_dict: all([feature[key] == value for key, value in search_dict.items()])
        features = []
        for feature in self.features:
            if search_fn(feature, kwargs):
                features.append(feature)
        return [self.seq_from_feature(f) for f in features]

    def reverse_complement(self):
        return self._seqrecord.reverse_complement()

    def seq_from_feature(self, feature):
        return Sequence(id=self._feature_id(feature), name=self._feature_name(feature), seq=feature.seq(self))

    def has_property(self, prop):
        return bool([x for x in self.properties if x.property == prop.property and x._type == prop._type])


class Contig(Sequence):
    meta = {'allow_inheritance': True,
            'index_cls': False,
            'collection': "contig_collection",
            'indexes': [
                'organism',
                {"fields": ["organism", "features.locus_tag"]},
                {"fields": ["features.identifier"]}
            ]}
    organism = StringField()
    organelle = StringField(required=False)
    seq_collection_id = ReferenceField(SeqCollection)
    bigseq = BinaryField()

    def __init__(self, **kwargs):
        '''
        '''
        super(Sequence, self).__init__(**kwargs)
        self._seq_init()
        self.size.unit = "bp"

    def gene(self, name):
        # gene_filter = lambda  x : all([ x[key] == value for key,value in params.items() ])
        for f in self.features:
            if f.has_alias(name):
                return f
        raise NotFoundException("Gene not found: " + name)

    # def seq_from_feature(self, whole_gene_feature, cds_feature=None, regulatory_feature=None):
    #     if not cds_feature:
    #         cds_feature = whole_gene_feature
    #     cds_feature = copy.copy(cds_feature)
    #     cds_feature.location = cds_feature.location.relative_to(whole_gene_feature.location)
    #     protein_coding_gene = ProteinCodingGene(id=self._feature_id(whole_gene_feature),
    #                                             name=self._feature_name(whole_gene_feature),
    #                                             seq=whole_gene_feature.seq(self))
    #     protein_coding_gene.add_cds(whole_gene_feature.identifier,
    #                                 cds_feature.location.start, cds_feature.location.end, cds_feature.location.strand,
    #                                 "")
    #     return protein_coding_gene
