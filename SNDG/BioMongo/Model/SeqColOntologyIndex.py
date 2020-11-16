from mongoengine.document import Document
from mongoengine.fields import StringField, ListField, ReferenceField, IntField

from SNDG.BioMongo.Model.SeqCollection import SeqCollection


class SeqColOntologyIndex(Document):
    meta = {'allow_inheritance': True,
            'collection': "col_ont_idx",
            'indexes': [
                'keywords',
                {"fields": ["seq_collection_name", "term"]},
                {"fields": ["seq_collection_id", "term"]},
                {"fields": ["seq_collection_name", "keywords"]},
                {"fields": ["seq_collection_id", "keywords"]}
            ]}

    term = StringField(required=True)
    name = StringField()
    count = IntField(required=True, default=0)
    order = IntField(default=0)
    keywords = ListField(StringField(), default=[])
    ontology = StringField(required=True)
    database = StringField(required=False)
    seq_collection_name = StringField(required=True)
    seq_collection_id = ReferenceField(SeqCollection)


class VarSeqColOntologyIndex(Document):
    meta = {'allow_inheritance': True,
            'collection': "var_col_ont_idx",
            'indexes': [
                'keywords',
                {"fields": ["seq_collection_name", "term"]},
                {"fields": ["seq_collection_id", "term"]},
                {"fields": ["seq_collection_name", "keywords"]},
                {"fields": ["seq_collection_id", "keywords"]}
            ]}

    term = StringField(required=True)
    name = StringField()
    count = IntField(required=True, default=0)
    order = IntField(default=0)
    keywords = ListField(StringField(), default=[])
    ontology = StringField(required=True)
    database = StringField(required=False)
    seq_collection_name = StringField(required=True)
    seq_collection_id = ReferenceField(SeqCollection)

#     @queryset_manager
#     def objects(doc_cls, queryset):
#         # This may actually also be done by defining a default ordering for
#         # the document, but this illustrates the use of manager methods
#         return queryset
