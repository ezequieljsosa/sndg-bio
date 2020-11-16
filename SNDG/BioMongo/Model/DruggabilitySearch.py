from mongoengine import DynamicEmbeddedDocument
from mongoengine.fields import StringField, ListField, EmbeddedDocumentField


class DruggabilitySearch(DynamicEmbeddedDocument):
    meta = {'allow_inheritance': True, 'abstract': False, 'strict': False}
    name = StringField()

    def __init__(self, **kwargs):
        '''
        '''
        super(DynamicEmbeddedDocument, self).__init__(**kwargs)


class PocketDruggabilitySearch(DruggabilitySearch):
    pocket = StringField(required=True)


class StructureDruggabilitySearch(DruggabilitySearch):
    structure = StringField(required=True)
    pockets = ListField(EmbeddedDocumentField(PocketDruggabilitySearch), default=[])


class ProteinDruggabilitySearch(DruggabilitySearch):
    meta = {'allow_inheritance': True, 'abstract': False, 'strict': False}
    structures = ListField(EmbeddedDocumentField(StructureDruggabilitySearch), default=[])
