
from datetime import datetime
from peewee import Proxy, Model, IntegerField, ForeignKeyField, CharField, BooleanField, DateTimeField, DoubleField, \
    AutoField, TextField


sqldb = Proxy()

class Pangenome(Model):
    id = AutoField(primary_key=True)
    name = CharField()
    kingdom = CharField()
    organelles = IntegerField(null=True)
    plasmids = IntegerField(null=True)
    chromosomes = IntegerField(null=True)
    tax = CharField(null=True)
    deprecated = BooleanField(default=False)
    date = DateTimeField(default=datetime.now)
    ncbi_id = IntegerField(null=True)
    class Meta:
        database = sqldb

class Strain(Model):
    id = AutoField(primary_key=True)
    pangenome = ForeignKeyField(Pangenome, related_name='genomes',
                          db_column="pangenome_id")
    name = CharField(null=True)
    acc = CharField()
    tax = CharField(null=True)
    status = CharField(null=True)
    deprecated = BooleanField(default=False)
    date = DateTimeField(default=datetime.now)
    loaded = BooleanField()
    ncbi_id = IntegerField(null=True)

    class Meta:
        database = sqldb