'''
Created on Oct 13, 2015

@author: eze
'''
from mongoengine.document import DynamicDocument
from mongoengine.fields import StringField, ListField


class Ontology(DynamicDocument):
    '''
    classdocs
    '''
    meta = {'allow_inheritance': True,
            'collection': "ontologies",
            'strict': False,
            'indexes': [
                'term', 'keywords',
                {"fields": ["ontology", "term"]},
                {"fields": ["ontology", "keywords"]}
            ]}

    term = StringField(required=True)
    name = StringField()
    ontology = StringField(required=True)
    databases = ListField(StringField(), default=[])
    database = StringField(default="")
    description = StringField(default="")
    keywords = ListField(StringField(), default=[])
    parents = ListField(StringField(), default=[])

    children = ListField(StringField(), default=[])
    successors = ListField(StringField(), default=[])
    subclases = ListField(StringField(), default=[])

    successors_relationships = ListField(StringField(), default=[])

    def __repr__(self):
        return "Ontology(" + ",".join([str(k) + "=" + v.__repr__() for k, v in self._data.items()]) + ")"

    def __str__(self):
        return "Ontology(" + ",".join([str(k) + "=" + str(v) for k, v in self._data.items()]) + ")"
