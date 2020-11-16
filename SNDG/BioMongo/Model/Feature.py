from mongoengine import EmbeddedDocument
from mongoengine.fields import StringField, ListField, EmbeddedDocumentField, DynamicField, ObjectIdField, IntField, \
    DictField

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq

from SNDG.BioMongo.Model.Alignment import SimpleAlignment


class Location(EmbeddedDocument):
    meta = {'allow_inheritance': True}
    base = StringField(max_length=50, required=False)
    start = IntField(required=True)
    end = IntField(required=True)
    strand = IntField(required=False)

    def range(self):
        return range(self.start, self.end)

    def relative_to(self, otherLocation):
        return self.__class__(base=self.base, start=self.start - otherLocation.start,
                              end=self.end - otherLocation.start)

    def __len__(self):
        return abs(self.end - self.start)

    def __str__(self):
        return self.base + ":" if self.base else "" + str(self.start) + "-" + str(self.end)


class DNALocation(Location):
    strand = IntField(required=True)

    def relative_to(self, otherLocation):
        loc = Location.relative_to(self, otherLocation)
        loc.strand = self.strand
        return loc

    def sense(self, sequence):
        if self.strand == 1:
            return sequence.seq[self.start:self.end]
        else:
            return str(Seq(sequence.seq[self.start:self.end]).reverse_complement())

    def __str__(self):
        return self.base + ":" if self.base else "" + str(self.start) + "-" + str(self.end) + "(" + (
            "+" if self.strand == 1 else "-") + ")"


class Feature(EmbeddedDocument):
    _id = ObjectIdField()  # TODO poner required=True
    source = StringField(required=False)
    evidence = StringField()
    identifier = StringField()
    location = EmbeddedDocumentField(Location)
    type = StringField(max_length=30, required=True)
    features = ListField(DynamicField())
    locus_tag = StringField(required=False)
    alias = ListField(StringField(), default=[])
    aln = EmbeddedDocumentField(SimpleAlignment)
    qualifiers = DictField(required=False)

    meta = {'allow_inheritance': True}

    def __init__(self, **kwargs):
        '''
        '''

        super(EmbeddedDocument, self).__init__(**kwargs)
        self._seq_feature = SeqFeature(FeatureLocation(self.location.start, self.location.end),
                                       ref=self.location.base, type=self.type, ref_db=self.identifier)
        self._instance = None
        self.features = []

    def seq(self, sequence):
        return sequence.seq[self.location.start:self.location.end]

    def _compound_name(self, other):
        return self.identifier + "_" + other.identifier

    def __contains__(self, other):
        if self.location.base == other.location.base:
            return self.location.start <= other.location.start and self.location.end >= other.location.end
        return False

    def __and__(self, other):
        if self.location.base == other.location.base:
            intersect = sorted(list(set(range(self.location.start, self.location.end)) & set(
                range(other.location.start, other.location.end))))
            if intersect:
                return Feature(location=Location(start=intersect[0], end=intersect[-1]),
                               source="operation", type="intersect",
                               identifier=self._compound_name(other))
            else:
                return Feature(location=Location(start=-1, end=-1),
                               source="operation", type="intersect",
                               identifier=self._compound_name(other))
        else:
            raise Exception("bases do not match")

    def __str__(self):
        return "Feature(ref={ref},location={location}, features={features_count}) ".format(
            location=str(self.location), features_count=str(len(self.features)), ref=self.location.base)

    def has_alias(self, name):
        if self.identifier == name:
            return True
        if hasattr(self, "alias"):
            return name in self.alias
        return False

    def __len__(self):
        return len(self.location)
