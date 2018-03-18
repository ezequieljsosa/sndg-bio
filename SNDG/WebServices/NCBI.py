'''
Created on Feb 21, 2017

@author: eze
'''

import datetime
import logging
import os
import subprocess as sp
import sys

from peewee import ForeignKeyField, CharField, Model, BooleanField, \
    DateTimeField, CompositeKey
from peewee import MySQLDatabase
from tqdm import tqdm

from Bio import Entrez
from SNDG import execute, init_log

mysql_db = MySQLDatabase('sndg', user="root", password="mito")

ftp_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all"


class Submitter(Model):
    '''
    classdocs
    '''

    # parent = ForeignKeyField('self', related_name='children')
    name = CharField()
    source = CharField()  # ncbi embl bold ...
    rejected = BooleanField(default=False)
    created = DateTimeField(default=datetime.datetime.now)
    modified = DateTimeField()

    def save(self, *args, **kwargs):
        self.modified = datetime.datetime.now()
        return super(Submitter, self).save(*args, **kwargs)

    class Meta:
        database = mysql_db

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return str(self.__data__)


class ExternalResource(Model):
    '''
    classdocs
    '''
    # id = IntegerField(primary_key=True)
    # submitter = ForeignKeyField(Submitter)
    name = CharField()
    identifier = CharField()
    type = CharField()  # ["gene","genome","bioproject","biosample","dbvar","sra","pubmed","protein","nuccore","assembly"]
    created = DateTimeField(default=datetime.datetime.now)
    modified = DateTimeField()

    def save(self, *args, **kwargs):
        self.modified = datetime.datetime.now()
        return super(ExternalResource, self).save(*args, **kwargs)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return str(self.__data__)

    class Meta:
        database = mysql_db
        indexes = (
            (('submitter', 'identifier',), True),
        )


class ExternalAssembly(ExternalResource):
    '''
    Se puede relacionar con sus bioentry haciendo:
    SELECT be.bioentry_id from bioentry be, bioentry_dbxref xr, dbxref r
    WHERE r.accession = "GCF_001742375.1"
        AND r.dbxref_id = xr.dbxref_id
        AND be.bioentry_id = xr.bioentry_id;

    '''
    assembly_accession = CharField()
    assembly_name = CharField()
    genome = CharField()

    class Meta:
        database = mysql_db
        indexes = (
            (('submitter', 'identifier',), True),
        )

    def ftpurl(self, dtype):
        return "/".join([ftp_url
                            , self.assembly_accession[0:3], self.assembly_accession[4:7]
                            , self.assembly_accession[7:10], self.assembly_accession[10:13]
                            , self.assembly_accession + "_" + self.assembly_name.replace(" ", "_").replace("#", "_")
                            , self.assembly_accession + "_" + self.assembly_name.replace(" ", "_").replace("#",
                                                                                                           "_") + "_" + dtype
                         ])

    def file_name(self, dtype):
        return self.assembly_accession + "_" + self.assembly_name.replace(" ", "_").replace("#", "_") + "_" + dtype

    def gbk_file_name(self):
        return self.file_name("genomic.gbff.gz")

    def gff_file_name(self):
        return self.file_name("genomic.gff.gz")

    def genome_file_name(self):
        return self.file_name("genomic.fna.gz")

    def prots_file_name(self):
        return self.file_name("protein.faa.gz")

    def dowload_from_ftp(self, dtype, dst_dir):
        current = os.getcwd()
        if not os.path.exists(dst_dir):
            os.makedirs(dst_dir)
        os.chdir(dst_dir)
        cmd = "ftp_proxy=http://proxy.fcen.uba.ar:8080/ wget --timeout=20 -q " + self.ftpurl(dtype)
        sp.call(cmd, shell=True)
        sp.call("gunzip " + self.file_name(dtype), shell=True)
        os.chdir(current)
        return dst_dir + "/" + self.file_name(dtype).replace(".gz", "")

    def download_gff(self, dst_dir):
        return self.dowload_from_ftp("genomic.gff.gz", dst_dir)

    def download_gbk(self, dst_dir):
        return self.dowload_from_ftp("genomic.gbff.gz", dst_dir)

    def download_genome(self, dst_dir):
        return self.dowload_from_ftp("genomic.fna.gz", dst_dir)

    def download_prots(self, dst_dir):
        return self.dowload_from_ftp("protein.faa.gz", dst_dir)


class AssemblySubmitters(Model):
    submitter = ForeignKeyField(Submitter)
    resource = ForeignKeyField(ExternalAssembly)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return str(self.__data__)

    class Meta:
        database = mysql_db
        primary_key = CompositeKey('submitter', 'resource')


_log = logging.getLogger(__name__)


class NCBIAssembly():

    def name(self, data):
        return data["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyName"]

    def attributes(self, data):
        return data["DocumentSummarySet"]["DocumentSummary"][0]


class NCBIProject():

    def name(self, data):
        return data["DocumentSummarySet"]["DocumentSummary"][0]["Project_Name"]


class NCBIGene():

    def name(self, data):
        return data["DocumentSummarySet"]["DocumentSummary"][0]["NomenclatureName"]


class NCBIBiosample():

    def name(self, data):
        return data["DocumentSummarySet"]["DocumentSummary"][0]["Title"]


class NCBIReads():
    def name(self, data):
        import xmltodict
        return xmltodict.parse("<x>" + data[0]["ExpXml"] + "</x>")["x"]["Summary"]["Title"]


class NCBIPubmed():

    def name(self, data):
        return data[0]["Title"]


class NCBIProtein():

    def name(self, data):
        return data[0]["Title"]


class NCBINucleotide():

    def name(self, data):
        return data[0]["Title"]


class NCBIGenome():

    def name(self, data):
        return data[0]["DefLine"]


class AssemblyNotFoundError(Exception):

    def __init__(self, assembly_name):
        self.assembly_name = assembly_name


class NCBI(object):
    f_mRNA = "mRNA"
    f_CDS = "CDS"
    ftypes = ["rRNA", "ncRNA", f_mRNA, "gene", "exon", f_CDS, "rRNA", "tRNA", "tmRNA"]

    # https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly
    dbs = ["bioproject", "biosample", "dbvar", "sra", "assembly"]  # "nuccore","protein",
    # "genome"
    dbs_con_submitter = ["assembly"]  # "bioproject", "biosample",

    # "pubmed", ?
    # "gene", ?

    def __init__(self):
        '''
        Constructor
        '''

        self.resource_handler = {
            "assembly": NCBIAssembly(), "bioproject": NCBIProject(), "biosample": NCBIBiosample(),
            "gene": NCBIGene(), "genome": NCBIGenome(), "dbvar": None, "sra": NCBIReads(), "pubmed": NCBIPubmed(),
            "protein": NCBIProtein(), "nuccore": NCBINucleotide()

        }

    def obtener_proyectos_argentinos(self):
        all_arg_proj = Entrez.read(Entrez.esearch(db="bioproject", term="argentina", retmax=10000))
        _log.info("existen %s proyectos con la palabra argentina" % all_arg_proj["Count"])
        for pid in tqdm(all_arg_proj["IdList"]):
            yield Entrez.read(Entrez.esummary(db="bioproject", id=pid))["DocumentSummarySet"]["DocumentSummary"][0]

    def obtener_submitters(self):
        submitters = []
        for bioproject in self.obtener_proyectos_argentinos():
            submitters = submitters + bioproject["Submitter_Organization_List"]

        return set(submitters)

    def actualizar_submitters(self, submitters):

        for submitter_raw in tqdm(submitters):
            submitter = str(submitter_raw.encode("utf-8"))
            submitter_in_db = Submitter.select().where((Submitter.source == "ncbi") & (Submitter.name == submitter))
            if not submitter_in_db.exists():
                submitter_in_db = Submitter(name=submitter, source="ncbi")
                submitter_in_db.save()

    def download(self, accesion, db, dst, dstformat):
        cmd = 'wget -O %s "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=%s&report=%s&sort=&id=%s&from=begin&to=end"'
        params = (dst, db, dstformat, accesion)

        execute(cmd % params, shell=True)

    def update_assemblies(self, submitters):
        db = "assembly"
        for submitter in tqdm(submitters):
            data = Entrez.read(Entrez.esearch(db=db, term='"' + submitter.name + '"[Submitter Organization]'))
            for recurso_id in data["IdList"]:
                er_in_db = ExternalAssembly.select().where(ExternalAssembly.identifier == recurso_id)
                if not er_in_db.exists():
                    recurso = Entrez.read(Entrez.esummary(db=db, id=recurso_id))
                    data = self.resource_handler[db].attributes(recurso)
                    name = str(self.resource_handler[db].name(recurso))
                    genome = str(data["SpeciesName"])
                    ea = ExternalAssembly(type=db, name=name, identifier=recurso_id
                                          , assembly_accession=data['AssemblyAccession'], genome=genome
                                          , assembly_name=data['AssemblyName'])
                    ea.save(force_insert=True)
                else:
                    ea = list(er_in_db)[0]

                AssemblySubmitters.create(submitter=submitter, resource=ea)

    def actualizar_recursos_externos(self):
        submitters = Submitter.select().where((Submitter.source == "ncbi") & (Submitter.rejected == False))
        for submitter in submitters:
            for db in NCBI.dbs_con_submitter:
                data = Entrez.read(Entrez.esearch(db=db, term='"' + submitter.name + '"[Submitter Organization]'))
                if 'ErrorList' in data:
                    _log.warn("error buscando %s en %s : %s " % (submitter.name, db, str(data['ErrorList'])))
                if db == "bioproject":
                    bioprojects = data["IdList"]
                if db == "assembly":
                    assemblies = data["IdList"]
                for recurso_id in data["IdList"]:
                    er_in_db = ExternalResource.select().where((ExternalResource.submitter == submitter) &
                                                               (ExternalResource.identifier == recurso_id))
                    if not er_in_db.exists():
                        recurso = Entrez.read(Entrez.esummary(db=db, id=recurso_id))
                        nombre = str(self.resource_handler[db].name(recurso))
                        ExternalResource(submitter=submitter, type=db, name=nombre, identifier=recurso_id).save()
                        if db == "assembly":
                            ExternalAssembly(submitter=submitter, type=db, name=nombre, identifier=recurso_id
                                             , assembly_accession=recurso['Assembly_Accession']
                                             , assembly_name=recurso['Assembly_Name']).save()

            for db in (set(NCBI.dbs) - set(NCBI.dbs_con_submitter)):
                for bioproject in bioprojects:
                    data = Entrez.read(Entrez.esearch(db=db, term=bioproject + '[BioProject]'))
                    if 'ErrorList' in data and not data['ErrorList']['PhraseNotFound']:
                        _log.warn("error buscando %s en %s : %s " % (submitter, db, str(data['ErrorList'])))
                    for recurso_id in data["IdList"]:
                        er_in_db = ExternalResource.select().where((ExternalResource.submitter == submitter) &
                                                                   (ExternalResource.identifier == recurso_id))
                        if not er_in_db.exists():
                            recurso = Entrez.read(Entrez.esummary(db=db, id=recurso_id))
                            nombre = self.resource_handler[db].name(recurso).encode("utf-8")
                            ExternalResource(submitter=submitter, type=db, name=nombre, identifier=recurso_id).save()

            db = "genome"
            for assembly in assemblies:
                data = Entrez.read(Entrez.esearch(db=db, term=assembly + '[AssemblyID]'))
                if 'ErrorList' in data and not data['ErrorList']['PhraseNotFound']:
                    _log.warn("error buscando %s en %s : %s " % (submitter, db, str(data['ErrorList'])))
                for recurso_id in data["IdList"]:
                    er_in_db = ExternalResource.select().where((ExternalResource.submitter == submitter) &
                                                               (ExternalResource.identifier == recurso_id))
                    if not er_in_db.exists():
                        recurso = Entrez.read(Entrez.esummary(db=db, id=recurso_id))
                        nombre = self.resource_handler[db].name(recurso).encode("utf-8")
                        # 'Assembly_Name': 'v.1.0'
                        # 'Assembly_Accession': 'GCA_000004255.1'
                        ExternalResource(submitter=submitter, type=db, name=nombre, identifier=recurso_id).save()

    def assembly_name_from_acc(self, assembly_accession, dst="/tmp/assembly.gb"):
        # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term=GCA_000009645.1
        esearch = Entrez.read(Entrez.esearch(db="assembly", term=assembly_accession))

        if esearch["IdList"]:
            assembly_id = esearch["IdList"][0]
            #             https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id=31148
            summary = Entrez.esummary(db="assembly", id=assembly_id)
            summary = Entrez.read(summary, validate=False)
            return summary["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyName"]
        else:
            raise AssemblyNotFoundError(assembly_accession)


def org_name(org_complete_name):
    arr = org_complete_name.replace("sequence", "").replace("genome", "").replace(",", "").replace(".", "").replace(" ",
                                                                                                                    "").replace(
        ":", "").split(" ")
    if len(arr) == 1:
        name = org_complete_name
    elif len(arr) == 2:
        name = arr[0][0:4] + arr[1][0:4]
    else:
        name = arr[0][0:4] + arr[1][0:4]
        if "subsp" in arr:
            if "strain" in arr:
                name += "".join(arr[arr.index("subsp") + 1: arr.index("strain")][:2])
            else:
                name += "".join(arr[arr.index("subsp") + 1: arr.index("subsp") + 3])
        if "strain" in arr:
            name += "".join(arr[arr.index("strain") + 1: arr.index("strain") + 3])
        else:
            name += arr[-1] if arr[-1] not in arr[0:2] else ""
    return name


if __name__ == "__main__":
    logger = logging.getLogger('peewee')
    logger.setLevel(logging.INFO)
    init_log()
    # Submitter.create_table()
    # ExternalResource.create_table()
    ExternalAssembly.create_table()
    AssemblySubmitters.create_table()

    Entrez.email = "ezejajaja@hotmail.com"
    # NCBI().actualizar_submitters(NCBI().obtener_submitters())

    submitters = Submitter.select().where((Submitter.source == "ncbi") & (Submitter.rejected == False))
    NCBI().update_assemblies(submitters)

    sys.exit()
    nombres = []
    # s = Submitter.select(Submitter.name == "sndg").get()

    # usehistory="y"
    buscar = True
    init = 0
    while buscar:
        retstart = init
        pepe = Entrez.read(Entrez.esearch(db="assembly", retmax="200",
                                          term='"representatives"[Filter]AND"latest refseq"[Filter]AND"bacteria"[Filter]AND"complete genome"[Assembly Level]'
                                          , usehistory="Y", retstart=str(retstart)))
        if pepe["IdList"]:
            juan = Entrez.read(Entrez.esummary(db="assembly", id=",".join(pepe["IdList"])), validate=False)
            for assembly in juan["DocumentSummarySet"]["DocumentSummary"]:
                # genome = str(Entrez.read(Entrez.elink(dbfrom="assembly", id=assembly.attributes["uid"], linkname="assembly_genome")) [0]['IdList'][0])
                genome = str(assembly["SpeciesName"])
                assert assembly['LastMajorReleaseAccession'].startswith("GCF")
                if not ExternalAssembly.select().where(
                        ExternalAssembly.assembly_accession == str(assembly["AssemblyAccession"])).count():
                    ExternalAssembly(assembly_accession=str(assembly["AssemblyAccession"]),
                                     assembly_name=str(assembly["AssemblyName"]),
                                     genome=genome,
                                     submitter_id=74,
                                     name=genome + " ref assembly",
                                     type="assembly",
                                     identifier=str(assembly.attributes["uid"])
                                     ).save()
                    # "representatives"[Filter]  "latest refseq"[Filter] "bacteria"[Filter] "complete genome"[Assembly Level] "chromosome"[Assembly Level]
            init = init + 200
        else:
            buscar = False

#         workdir = "/data/projects/ncbi_dump/" 
#         os.chdir(workdir)
#         for x in ExternalAssembly.select():
#             aworkdir= str(x.id)  
#             file_name =  x.gbk_file_name().replace(".gz","")                        
#             if not os.path.exists(aworkdir):                
#                 os.makedirs( aworkdir )
#                 x.download_gbk(aworkdir)                
#             result = aworkdir + "/" + file_name
#             if os.path.exists(result):
#                 org =  org_name(list(bpio.parse(result,"gb"))[0].description)
#                 print org
#                 nombres.append(org)
#             else:
#                 print "no exists: " + x.identifier
#         print len(nombres)
#         print len(set(nombres))

from SNDG.BioMongo.Process.Taxon import Tax, tax_db
tax_db.initialize(MySQLDatabase('bioseqdb', user='root', passwd="mito"))
tax = Tax.getTax(1872703)
for col in db.sequence_collection.find({"tax":{"$exists":0}, "ncbi_assembly":{"$exists":1} },{"name":1,"ncbi_assembly":1,"_id":0}):
    esearch = Entrez.read(Entrez.esearch(db="assembly", term= col["ncbi_assembly"] ))["IdList"][0]
    summary = Entrez.esummary(db="assembly", id=esearch)
    summary = Entrez.read(summary, validate=False)
    juan = dict(summary["DocumentSummarySet"]["DocumentSummary"][0])
    #pepe["AssemblyStatus"]
    tax = Tax.getTax(str(juan["Taxid"]))
    tax = {
        "tid" : float(tax.ncbi_taxon_id),
                         "superkingdom" : [y for y in [x for x in Tax.parents(tax) if x.node_rank == "superkingdom"][0].names
                                           if y.name_class == "scientific name"  ][0].name  ,
        "name" : [x for x in tax.names if x.name_class == "scientific name"][0].name
        }
    print tax
    print db.sequence_collection.update({"name":col["name"]},{"$set":{"tax":tax}})
