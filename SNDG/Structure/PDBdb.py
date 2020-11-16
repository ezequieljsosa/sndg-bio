"""

"""
from datetime import datetime
from peewee import Proxy, Model, IntegerField, ForeignKeyField, CharField, BooleanField, DateTimeField, DoubleField, \
    AutoField, TextField

sqldb = Proxy()


class PDB(Model):
    id = AutoField(primary_key=True)
    code = CharField()
    resolution = DoubleField(default=20)
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
    resid = IntegerField()
    icode = CharField()

    type = CharField()  # HETATOM / ATOM
    disordered = BooleanField()

    has_ca = BooleanField(null=True, default=None)
    only_ca = BooleanField(null=True, default=None)

    class Meta:
        database = sqldb


class Atom(Model):
    id = AutoField(primary_key=True)

    serial = IntegerField()
    name = CharField()
    residue = ForeignKeyField(Residue, related_name='atoms',
                              db_column="residue_id")
    x = DoubleField()
    y = DoubleField()
    z = DoubleField()
    occupancy = DoubleField()
    bfactor = DoubleField()
    anisou = DoubleField(null=True)
    element = CharField()

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


class Property(Model):
    id = AutoField(primary_key=True)
    name = CharField()
    description = TextField(null=True)

    class Meta:
        database = sqldb


class PDBProperty(Model):
    id = AutoField(primary_key=True)
    pdb = ForeignKeyField(PDB, related_name='properties',
                          db_column="pdb_id")
    property = ForeignKeyField(PDB, db_column="property_id")

    value = DoubleField(null=True)
    tag = CharField(null=True)

    class Meta:
        database = sqldb


class ResidueProperty(Model):
    id = AutoField(primary_key=True)
    residue = ForeignKeyField(Residue, related_name='properties',
                              db_column="residue_id")
    property = ForeignKeyField(PDB, db_column="property_id")
    value = DoubleField(null=True)
    tag = CharField(null=True)

    class Meta:
        database = sqldb


class ResidueSetProperty(Model):
    id = AutoField(primary_key=True)
    residue_set = ForeignKeyField(ResidueSet, related_name='properties',
                                  db_column="residue_set_id")
    property = ForeignKeyField(PDB, db_column="property_id")
    value = DoubleField(null=True)
    tag = CharField(null=True)

    class Meta:
        database = sqldb


class ChainProperty(Model):
    id = AutoField(primary_key=True)
    pdb = ForeignKeyField(PDB, related_name='chain_props',
                          db_column="pdb_id")
    chain = CharField()
    property = ForeignKeyField(PDB, db_column="property_id")
    value = DoubleField(null=True)
    tag = CharField(null=True)

    class Meta:
        database = sqldb


class AtomProperty(Model):
    id = AutoField(primary_key=True)
    atom = ForeignKeyField(PDB, related_name='properties',
                           db_column="atom_id")
    property = ForeignKeyField(PDB, db_column="property_id")
    value = DoubleField(null=True)
    tag = CharField(null=True)

    class Meta:
        database = sqldb


if __name__ == "__main__":
    from peewee import MySQLDatabase

    mysql_db = MySQLDatabase('pdbdb', user="root", password="mito")
    sqldb.initialize(mysql_db)
    tables = [PDB, Residue, Atom, ResidueSet, ResidueSetResidue, PDBProperty,
              ResidueProperty, ResidueSetProperty, ChainProperty, AtomProperty]
    # for x in reversed(tables):
    #     x.drop_table()
    # for x in tables:
    #     x.create_table()

    # Property.create_table()
    # qmean_s = ['qr_solvation', 'qr_ss_agreement', 'qr_QMEANDisCo', 'qr_exposed', 'qr_cbeta', 'qr_all_atom', 'qr_QMEAN',
    #            'qr_acc_agreement', 'qr_dist_const', 'qr_torsion']
    # qmean_r = ['q_torsion_zscore', 'q_QMEAN6_zscore', 'q_cbeta_zscore', 'q_interaction_norm', 'q_acc_agreement_zscore',
    #            'q_QMEAN4_norm', 'q_torsion_norm', 'q_packing_zscore', 'q_interaction_zscore', 'q_QMEAN6_norm',
    #            'q_QMEAN4_zscore', 'q_residues', 'q_packing_norm', 'q_ss_agreement_norm', 'q_ss_agreement_zscore',
    #            'q_cbeta_norm', 'q_acc_agreement_norm']
    # fpocket = ['f_Volume score', 'f_Polarity score', 'f_Proportion of polar atoms', 'f_Druggability Score',
    #            'f_Cent of mass - Alpha Sphere max dist', 'f_Flexibility', 'f_Alpha sphere density', 'f_Charge score',
    #            'f_Apolar alpha sphere proportion', 'f_Number of Alpha Spheres', 'f_Volume', 'f_Hydrophobicity score',
    #            'f_Apolar SASA', 'f_Score', 'f_Mean local hydrophobic density', 'f_Polar SASA', 'f_Total SASA',
    #            'f_Mean alpha sphere radius', 'f_Mean alp sph solvent access']
    # for f in qmean_r + qmean_s + fpocket:
    #     Property(name=f).save()
