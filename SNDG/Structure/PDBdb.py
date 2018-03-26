"""

"""
from datetime import datetime
from peewee import Proxy, Model, IntegerField, ForeignKeyField, CharField, BooleanField, DateTimeField, DoubleField, \
    AutoField

sqldb = Proxy()


class PDB(Model):
    id = AutoField(primary_key=True)
    code = CharField()
    resolution = IntegerField(default=20)
    experiment = CharField(null=True)
    tax = CharField(null=True)
    deprecated = BooleanField(default=False)
    date = DateTimeField(default=datetime.now)

    class Meta:
        database = sqldb


class Residue(Model):
    id = AutoField(primary_key=True)
    pdb = ForeignKeyField(PDB, related_name='residues',
                          db_column="pdb_id")
    chain = CharField()
    resname = CharField()
    resid = CharField()
    icode = CharField()

    # x = DoubleField()
    # y = DoubleField()
    # z = DoubleField()
    # occupancy = DoubleField()
    # tempFactor = DoubleField()
    # element = CharField()
    # charge = CharField()

    type = CharField()  # HETATOM / ATOM
    has_ca = BooleanField(null=True, default=None)
    only_ca = BooleanField(null=True, default=None)

    class Meta:
        database = sqldb


class ResidueSet(Model):
    id = AutoField(primary_key=True)
    pdb = ForeignKeyField(PDB, related_name='residue_sets',
                          db_column="pdb_id")
    name = CharField()  # Pocket / CSA

    class Meta:
        database = sqldb


class ResidueSetResidue(Model):
    id = AutoField(primary_key=True)
    residue = ForeignKeyField(Residue, related_name='residue_sets',
                              db_column="residue_id")
    residue_set = ForeignKeyField(ResidueSet, related_name='residues',
                                  db_column="residue_set_id")
    name = CharField()  # Pocket / CSA

    class Meta:
        database = sqldb


class PDBProperty(Model):
    id = AutoField(primary_key=True)
    pdb = ForeignKeyField(PDB, related_name='properties',
                          db_column="pdb_id")
    name = CharField()  # sasa / qmean
    value = DoubleField(null=True)
    tag = CharField(null=True)

    class Meta:
        database = sqldb


class ResidueProperty(Model):
    id = AutoField(primary_key=True)
    residue = ForeignKeyField(Residue, related_name='properties',
                              db_column="residue_id")
    name = CharField()  # sasa / qmean
    value = DoubleField(null=True)
    tag = CharField(null=True)

    class Meta:
        database = sqldb


class ResidueSetProperty(Model):
    id = AutoField(primary_key=True)
    residue_set = ForeignKeyField(ResidueSet, related_name='properties',
                                  db_column="residue_set_id")
    name = CharField()  # druggability / qmean
    value = DoubleField(null=True)
    tag = CharField(null=True)

    class Meta:
        database = sqldb


class ChainProperty(Model):
    id = AutoField(primary_key=True)
    pdb = ForeignKeyField(PDB, related_name='chain_props',
                          db_column="pdb_id")
    chain = CharField()
    name = CharField()  # druggability / qmean
    value = DoubleField(null=True)
    tag = CharField(null=True)

    class Meta:
        database = sqldb


if __name__ == "__main__":
    from peewee import MySQLDatabase

    mysql_db = MySQLDatabase('pdbdb', user="root", password="mito")
    sqldb.initialize(mysql_db)
    for x in reversed([PDB, Residue, ResidueSet, ResidueSetResidue, PDBProperty,
              ResidueProperty, ResidueSetProperty, ChainProperty]):
        x.drop_table()
    for x in [PDB, Residue, ResidueSet, ResidueSetResidue, PDBProperty,
              ResidueProperty, ResidueSetProperty, ChainProperty]:
        x.create_table()
