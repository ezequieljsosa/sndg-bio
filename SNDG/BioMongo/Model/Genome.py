from mongoengine import EmbeddedDocument
from mongoengine.fields import EmbeddedDocumentField, ListField, DynamicField

from SNDG.BioMongo.Model import BioProperties
from SNDG.BioMongo.Model.SeqCollection import SeqCollection
from SNDG.BioMongo.Model.exceptions import NotFoundException


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
