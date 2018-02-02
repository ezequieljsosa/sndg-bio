'''
'''
from mongoengine.document import DynamicEmbeddedDocument,EmbeddedDocument
from mongoengine.fields import StringField, ListField, EmbeddedDocumentField


class BioProperties(DynamicEmbeddedDocument):
    '''
    '''

    def __str__(self):
        result = "DynDoc("
        for k, v in self._data.items():
            result = result + str(k) + "=" + str(v) + ","
        return result[0:-1] + ")"

    def __repr__(self, *args, **kwargs):
        return self.__str__()


class Cluster(EmbeddedDocument):
    name = StringField(max_length=50, required=True)
    type = StringField(max_length=50, required=True)
    parts = ListField(EmbeddedDocumentField(BioProperties), default=[])

    def __repr__(self):
        return "Cluster(" + self.name + ", size: " + str(len(self.parts)) + ")"

    def __str__(self):
        return self.__repr__()

    def __nonzero__(self):
        return len(self._data["parts"]) > 0
