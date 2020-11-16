"""

"""
from mongoengine import Document
from mongoengine.fields import StringField, ListField, BooleanField, IntField


class Taxonomy(Document):
    '''

    '''
    meta = {'allow_inheritance': True, 'strict': False, 'collection': "tax"}

    name = StringField(required=True)
    names = ListField(StringField(), default=[])
    parents = ListField(StringField(), default=[])
    keywords = ListField(StringField(), default=[])
    ncbi_taxon_id = IntField(required=True)
    rank = StringField(required=True)
    genetic_code = IntField(required=False)
    mito_genetic_code = IntField(required=False)
