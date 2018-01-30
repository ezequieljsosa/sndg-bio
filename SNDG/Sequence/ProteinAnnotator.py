from tqdm import tqdm

from peewee import MySQLDatabase, Model, IntegrityError

from peewee import Proxy
from peewee import CharField, BooleanField, DateTimeField, ForeignKeyField
from playhouse.shortcuts import OperationalError

from datetime import date, datetime
import json

from SNDG import grouper


def json_serial(obj):
    """JSON serializer for objects not serializable by default json code"""

    if isinstance(obj, (datetime, date)):
        serial = obj.isoformat()
        return serial
    raise TypeError("Type %s not serializable" % type(obj))


class PABase(Model):
    sqldb = Proxy()

    def __str__(self):
        return json.dumps(self._data, default=json_serial)

    def __repr__(self):
        return self.__str__()


class Uniprot(PABase):
    uniprot = CharField(primary_key=True)
    uniref90 = CharField(null=True)
    deprecated = BooleanField(default=False)
    timestamp = DateTimeField(default=datetime.now)
    '''
    
    '''

    class Meta:
        database = PABase.sqldb


class Mapping(PABase):
    uniprot = ForeignKeyField(Uniprot, related_name='mappings',
                              db_column="uniprot_fk")
    db = CharField()
    value = CharField()

    class Meta:
        database = PABase.sqldb


class MappingProperty(PABase):
    mapping = ForeignKeyField(Mapping, related_name='properties',
                              db_column="mapping_fk")
    db = CharField()
    value = CharField()

    class Meta:
        database = PABase.sqldb


class ProteinAnnotator:
    tables = [Uniprot, Mapping, MappingProperty]

    def __init__(self):
        pass

    def connect_to_db(self, database='unipmap', user='root', password='', engine=MySQLDatabase):
        """class MyRetryDB(OperationalError):
            def __init__(self,engine,**kwargs):
                super().__init__(engine,kwargs)
        """

        # mysqldb = MyRetryDB(database, user=user, password=password)
        mysqldb = MySQLDatabase(database, user=user, password=password)
        PABase.sqldb.initialize(mysqldb)

    def prepare_uniprot_mapping(self):
        with open("idmapping_filtered.dat", "w") as h1:
            with open("idmapping.dat") as h2:
                for l in tqdm(h2, total=1519700351):
                    tupla = l.split()
                if tupla[1] not in ["GI", "UniRef50", "UniRef100",
                                    "EMBL", "RefSeq", "RefSeq_NT", "CRC64", "GeneID", "UniProtKB-ID",
                                    "EMBL-CDS", "UniParc", "OMA", "EnsemblGenome", "STRING", "UniPathway",
                                    "Ensembl_TRS", "Ensembl_PRO", "Allergome", "PATRIC"]:
                    h1.write(l)

    def create_db(self):
        # PABase.sqldb.execute_sql('DROP DATABASE ' + PABase.sqldb.database + ";")
        # PABase.sqldb.execute_sql('CREATE DATABASE ' + PABase.sqldb.database + ";")

        for t in ProteinAnnotator.tables:
            t.create_table()

    def populate_sql(self, unip_id_mapping, goa_id_mapping):
        import subprocess as sp

        total = int(sp.check_output("wc -l " + unip_id_mapping, shell=True).split()[0])
        with open(unip_id_mapping) as h:
            all_unips = {}
            previous = None
            for lines in grouper(tqdm(h, total=total), 10000):
                data = []
                for line in lines:
                    if line:
                        unip, _, _ = line.strip().split("\t")
                        if unip != previous and unip not in all_unips:
                            data.append((unip,))
                            previous = unip
                        all_unips[unip] = 1

                fields = [Uniprot.uniprot]
                with PABase.sqldb.atomic():
                    op = Uniprot.insert_many(data, fields=fields)
                    op.execute()

        with open(unip_id_mapping) as h:
            for lines in grouper(tqdm(h, total=total), 10000):
                data = []
                for line in lines:
                    if line:
                        unip, db, db_id = line.strip().split("\t")
                        data.append((unip, db, db_id))

                fields = [Mapping.uniprot, Mapping.db, Mapping.value]

                with PABase.sqldb.atomic():
                    op = Mapping.insert_many(data, fields=fields)
                    op.execute()

        total = int(sp.check_output("wc -l " + goa_id_mapping, shell=True).split()[0])
        with open(goa_id_mapping) as h:
            for lines in grouper(tqdm(h, total=total), 10000):
                data_m = []
                for line in lines:
                    if line:
                        if line.startswith("!") or line.startswith("gpa-version"):
                            continue
                        db, db_object_id, qualifiers, go_id, db_ref, eco_ev, _with, interacting_taxon, date, \
                        assigned_by = line.strip().split("\t")[:10]
                        if db != "UniProtKB":
                            continue
                        data_m.append((db_object_id, "GO", go_id))
                        if _with.startswith("EC:"):
                            data_m.append((db_object_id, "EC", _with))
                fields = [Mapping.uniprot, Mapping.db, Mapping.value]
                with PABase.sqldb.atomic():
                    op = Mapping.insert_many(data_m, fields=fields)
                    op.execute()


if __name__ == '__main__':
    pa = ProteinAnnotator()
    pa.connect_to_db(password="mito")
    pa.create_db()
    pa.populate_sql("/data/uniprot/idmapping_filtered.dat",
                    "/data/uniprot/goa/goa_uniprot_all.gpa")
