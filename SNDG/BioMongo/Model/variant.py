'''
Created on Aug 25, 2017

@author: eze
'''
from mongoengine.base.fields import ObjectIdField
from mongoengine.document import EmbeddedDocument, Document
from mongoengine.fields import StringField, ListField, IntField, \
    EmbeddedDocumentField,  DictField
from Bia.Model.sequence import ProteinDruggabilitySearch, Feature


class SampleAllele(EmbeddedDocument):
    
     
    sample = StringField( required=True)    
    annotations = DictField()
    
    
        
    

class Allele(EmbeddedDocument):
    
    _id = ObjectIdField() 
    samples = ListField(EmbeddedDocumentField(SampleAllele))
    alt = StringField( required=True)
    
    
        
    variant_type = ListField(StringField(),default=[])
    aa_ref = StringField( required=False)
    aa_alt = StringField( required=False)
    aa_pos = IntField( required=True)
    
    feature_ref = ObjectIdField()
    feature = EmbeddedDocumentField(Feature)

class Variant(Document):    
       
    
    organism = StringField( required=True)
    contig = StringField( required=True)
    pos = IntField( required=True)
    gene = StringField( required=False)
    prot_ref = ObjectIdField()
    ref = StringField( required=False)
    sample_alleles = ListField(EmbeddedDocumentField(Allele))
    ontologies = ListField(StringField())
    search = EmbeddedDocumentField(ProteinDruggabilitySearch,default=None)
    
    
    
    
    