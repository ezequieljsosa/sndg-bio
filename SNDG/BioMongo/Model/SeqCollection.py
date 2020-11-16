"""

"""
import datetime

from mongoengine import Document, EmbeddedDocument
from mongoengine.fields import StringField, ListField, EmbeddedDocumentField, BooleanField, ObjectIdField, IntField, \
    DateTimeField, DictField, FloatField,DynamicField

from SNDG.BioMongo.Model import BioProperties
from SNDG.BioMongo.Model.Pathway import PathwaySumary
from SNDG.BioMongo.Model.SeqColDruggabilityParam import SeqColDruggabilityParam
from SNDG.BioMongo.Model.exceptions import NotFoundException


class Metric(EmbeddedDocument):
    name = StringField(max_length=50, required=True)
    description = StringField(required=False)
    value = FloatField(required=False)
    values = ListField(FloatField(), required=False)
    labels = ListField(DynamicField(), required=False)

    def __str__(self):
        return "Metric: %s (%s) %s" % (self.name, self.description, str(self.value))


class AnnotationPipelineResult(EmbeddedDocument):
    name = StringField(required=True)
    timestamp = DateTimeField(required=False, default=datetime.datetime.now)
    version = IntField(default=0)
    inputs = DictField()
    results = DictField()


class Strain(EmbeddedDocument):
    name = StringField( required=True)
    description = StringField(default="")
    latitude = FloatField()
    longitude = FloatField()
    date = DateTimeField(default=datetime.datetime.now)
    country = StringField(default="")
    region = StringField(default="")
    properties = ListField(DictField(required=False),default=[])
    user = StringField(default="demo")


class StrainProject(EmbeddedDocument):
    _id = ObjectIdField()
    name = StringField( required=True)
    description = StringField(default="")
    date = DateTimeField(default=datetime.datetime.now)
    strains = ListField(StringField())
    trees = ListField(DictField())
    user = StringField(default="demo")

class StrainProp(EmbeddedDocument):
    name = StringField( required=True)
    description = StringField(default="")
    category = StringField(default="")
    options = ListField(StringField(),default=[])
    user = StringField(default="demo")


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
    kegg = ListField(EmbeddedDocumentField(PathwaySumary), default=[])

    ec_index = BooleanField()
    go_index = BooleanField()
    auth = StringField()
    version = IntField(default=0)
    pipelines = ListField(EmbeddedDocumentField(AnnotationPipelineResult), default=[])
    druggabilityParams = ListField(EmbeddedDocumentField(SeqColDruggabilityParam), default=[])

    strains = ListField(EmbeddedDocumentField(Strain), default=[])
    strainProjects = ListField(EmbeddedDocumentField(StrainProject), default=[])
    strainsProps = ListField(EmbeddedDocumentField(StrainProp), default=[])

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


class Genome(SeqCollection):
    '''

    '''
    meta = {'allow_inheritance': True, 'strict': False}
    assembly = EmbeddedDocumentField(BioProperties)

    def add_contig(self, contig):
        contig.seq_collection_id = self
        self._sequences.append(contig)

    def gene_contig(self, name):
        for contig in self._sequences:
            try:
                gene = contig.gene(name)
                if gene:
                    return gene, contig
            except NotFoundException as _:
                pass
        raise NotFoundException(name)

    def gene(self, name):
        gene, contig = self.gene_contig(name)
        return contig.seq_from_feature(gene)

    def gene_feature(self, name):
        gene, _ = self.gene_contig(name)
        return gene

    def proteins(self, name):
        return self.gene(name).proteins()

    def __str__(self):
        return "SeqCollection( {name}, contigs={contigs} )".format(name=self.name, contigs=str(len(self._sequences)))
