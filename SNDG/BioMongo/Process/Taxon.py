'''
Created on May 22, 2017

@author: eze
'''

from peewee import Proxy, DecimalField, CharField, Model, PrimaryKeyField, \
    MySQLDatabase, ForeignKeyField, CompositeKey

tax_db = Proxy()


class Tax(Model):
    '''
    classdocs
    '''

    scientific_name = "scientific name"

    taxon_id = PrimaryKeyField(db_column="taxon_id")
    ncbi_taxon_id = DecimalField()
    parent = ForeignKeyField('self', null=True, related_name="children", db_column="parent_taxon_id")
    node_rank = CharField()
    genetic_code = DecimalField()
    mito_genetic_code = DecimalField()
    left_value = DecimalField()
    right_value = DecimalField()

    def __str__(self):
        return "Tax(" + "".join([x + "=" + str(y) + "," for x, y in self.__data__.items()]) + ")"

    @classmethod
    def parents(cls, tax,return_self=True):
        if tax.taxon_id == 1:
            return []
        else:
            return ([tax] if return_self else []) + Tax.parents(tax.parent)


    @classmethod
    def getTax(cls, tax_id):
        return cls.select().where(Tax.ncbi_taxon_id == tax_id).get()

    class Meta:
        database = tax_db
        db_table = "taxon"


class TaxName(Model):
    '''
    classdocs
    '''



    taxon = ForeignKeyField(Tax, related_name="names", db_column="taxon_id")
    name = CharField(max_length=255)
    name_class = CharField(max_length=32)

    def __str__(self):
        return "TaxName(" + "".join([x + "=" + str(y) + "," for x, y in self.__data__.items()]) + ")"

    class Meta:
        database = tax_db
        db_table = "taxon_name"
        primary_key = CompositeKey('taxon_id', 'name', 'name_class')


if __name__ == '__main__':
    tax_db.initialize(MySQLDatabase('bioseqdb', user='root', passwd="mito"))
    tax = Tax.getTax(1872703)
    for x in tax.names:
        print x
    for t in Tax.parents(tax):
        print t
