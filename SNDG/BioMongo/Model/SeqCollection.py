"""

"""
import datetime

from mongoengine import Document, EmbeddedDocument
from mongoengine.fields import StringField, ListField, EmbeddedDocumentField, BooleanField, ObjectIdField, IntField, \
    DateTimeField, DictField, FloatField

from SNDG.BioMongo.Model.SeqColDruggabilityParam import SeqColDruggabilityParam
from SNDG.BioMongo.Model.Pathway import PathwaySumary


class Metric(EmbeddedDocument):
    name = StringField(max_length=50, required=True)
    description = StringField(required=False)
    value = FloatField(required=False)
    values = ListField(FloatField(), required=False)
    labels = ListField(StringField(), required=False)

    def __str__(self):
        return "Metric: %s (%s) %s" % (self.name, self.description, str(self.value))


class AnnotationPipelineResult(EmbeddedDocument):
    name = StringField(required=True)
    timestamp = DateTimeField(required=False, default=datetime.datetime.now)
    version = IntField(default=0)
    inputs = DictField()
    results = DictField()


class Strain(EmbeddedDocument):
    name = StringField(max_length=50, required=True)
    description = StringField(default="")
    latitude = FloatField()
    longitude = FloatField()
    date = DateTimeField(default=datetime.datetime.now)
    country = StringField(default="")
    region = StringField(default="")
    properties = DictField(required=False)
    projects = ListField(ObjectIdField(), default=[])


class TaxEDoc(EmbeddedDocument):
    tid = FloatField()
    superkingdom = StringField()
    name = StringField()


class DataUpload(EmbeddedDocument):
    name = StringField()
    timestamp = DateTimeField(default=datetime.datetime.now)
    uploader = StringField()
    errors = ListField(StringField())
    properties = ListField(StringField())


class SeqCollection(Document):
    '''

    '''
    meta = {'allow_inheritance': True, 'strict': False, 'collection': "sequence_collection"}

    name = StringField(max_length=50, required=True)
    type = StringField(required=False)
    ncbi_assembly = StringField(required=False)
    description = StringField(required=False, default="")
    organism = StringField(max_length=100)
    pathways = ListField(EmbeddedDocumentField(PathwaySumary), default=[])
    ec_index = BooleanField()
    go_index = BooleanField()
    auth = ObjectIdField()
    version = IntField(default=0)
    pipelines = ListField(EmbeddedDocumentField(AnnotationPipelineResult), default=[])
    druggabilityParams = ListField(EmbeddedDocumentField(SeqColDruggabilityParam), default=[])

    strainsProps = ListField(EmbeddedDocumentField(Strain), default=[])

    tax = EmbeddedDocumentField(TaxEDoc)
    statistics = ListField(EmbeddedDocumentField(Metric))
    uploads = ListField(EmbeddedDocumentField(DataUpload), default=[])

    def __init__(self, **kwargs):
        '''
        '''
        super(Document, self).__init__(**kwargs)
        self._sequences = []

    def druggabilityParam(self, name, uploader="demo"):
        return [x for x in self.druggabilityParams if x.name == name and x.uploader == uploader]

    def has_druggability_param(self, name, uploader="demo"):
        return bool(self.druggabilityParam(name, uploader))

    def add_drugability_props_to_genome(self, name, description, target, _type, options=None, user="demo"):
        """
        ptype: SeqColDruggabilityParamTypes
        """

        if not self.has_druggability_param(name):
            dp = SeqColDruggabilityParam(name=name, description=description, target=target,
                                         type=_type, uploader=user)
            if options:
                dp.options = options
            self.druggabilityParams.append(dp)
