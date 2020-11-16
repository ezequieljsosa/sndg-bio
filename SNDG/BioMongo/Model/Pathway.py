'''
Created on Sep 24, 2015

@author: eze
'''
from enum import Enum
from mongoengine.document import EmbeddedDocument
from mongoengine.fields import StringField, ListField, EmbeddedDocumentField, \
    IntField,DictField


class ChokepointType(Enum):
    consuming = "consuming"
    production = "production"
    double = "double"


class ChemSpecie(EmbeddedDocument):
    name = StringField( required=True)
    stoichiometry = IntField( )
    producers =  ListField(StringField())
    consumers =  ListField(StringField())
    

class ChockePoint(EmbeddedDocument):
    type = StringField( )
    specie = ChemSpecie( )
    

class Reaction(EmbeddedDocument):

    pathways = ListField(StringField( required=True))
    name = StringField( required=True)
    description = StringField( )
    substrates =  ListField(EmbeddedDocumentField(ChemSpecie))
    products = ListField(EmbeddedDocumentField(ChemSpecie))
    chockepoint = EmbeddedDocumentField(ChockePoint)
    
    def product(self,**kwargs):
        filtered = [x for x in self.products if all([ x[k] == v for k,v in kwargs.items()  ]) ]
        if filtered:
            return filtered[0]
        return None
    
    def substrate(self,**kwargs):
        filtered = [x for x in self.substrates if all([ x[k] == v for k,v in kwargs.items()  ]) ]
        if filtered:
            return filtered[0]
        return None

class PathwaySumary(EmbeddedDocument):
    term = StringField( required=True)
    name = StringField( required=True)
    count = IntField( required=True)
    properties = DictField(default={})
    
    def __str__(self):
        return "term='%s',name='%s',count='%i'" % (self.term, self.name,self.count)

    