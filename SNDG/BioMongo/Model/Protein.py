from mongoengine.fields import StringField, ListField, EmbeddedDocumentField, ObjectIdField, ReferenceField, FloatField, \
    DictField,EmbeddedDocument

from SNDG.BioMongo.Model.Sequence import Sequence
from SNDG.BioMongo.Model.SeqCollection import SeqCollection
from SNDG.BioMongo.Model.Pathway import Reaction
from SNDG.BioMongo.Model.Feature import Feature
from SNDG.Sequence.so import SO_TERMS


class ChEMBLAssay(EmbeddedDocument):
    assay = StringField(required=True)
    type = StringField(required=True)
    description = StringField(required=False)
    activities = ListField(DictField(), default=[])

class ChEMBL(EmbeddedDocument):
    target = StringField(required=True)
    assays = ListField(EmbeddedDocumentField(ChEMBLAssay), default=[])

class Protein(Sequence):
    meta = {'allow_inheritance': True, 'collection': "proteins",
            'index_cls': False,
            'indexes': [
                {"fields": ["organism", "name"]},
                {"fields": ["organism", "gene"]},
                {"fields": ["seq_collection_id", "keywords"]},
                {"fields": ["organism", "keywords"]},
                {"fields": ["organism", "alias"]},
                'keywords'
            ]}

    organism = StringField()
    gene = ListField(StringField(), required=True)
    gene_id = ObjectIdField()
    seq_collection_id = ReferenceField(SeqCollection)
    reactions = ListField(EmbeddedDocumentField(Reaction), default=[])
    reactome = ListField(DictField(), default=[])
    druggability = FloatField()
    tregulation = DictField(required=False)
    alias = ListField(StringField(), default=[])
    chembl = EmbeddedDocumentField(ChEMBL)

    def __init__(self, **kwargs):
        '''
        '''
        super(Sequence, self).__init__(**kwargs)
        self._seq_init()
        self.size.unit = "aa"

    def add_domain(self, name, start, end, source="unknown"):
        dn_feature = Feature(source=source, identifier=name, type=SO_TERMS["polypeptide_domain"],
                             location=Location(base=self.id, start=start, end=end))
        self.features.append(dn_feature)

    def domains(self):
        return [feature for feature in self.features if feature.type == SO_TERMS["polypeptide_domain"]]

    def domain_seqs(self):
        return self.seqs_from_features(type=SO_TERMS["polypeptide_domain"])

    def reaction(self, **kwargs):
        filtered_reactions = [x for x in self.reactions if all([x[k] == v for k, v in kwargs.items()])]
        if filtered_reactions:
            return filtered_reactions[0]
        return None

    def __str__(self):
        return self.name + " " + str(len(self.seq)) + "aa"

    def homologous(self):
        return [feature for feature in self.features if feature.type == SO_TERMS["homologous"]]
