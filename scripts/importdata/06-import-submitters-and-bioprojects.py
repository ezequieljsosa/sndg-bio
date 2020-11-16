import pymongo

from Bio import Entrez
from tqdm import tqdm
from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
from SNDG.WebServices.NCBI import ExternalAssembly, mysql_db, Submitter, BioProject, NCBI, AssemblySubmitters
from peewee import MySQLDatabase

mysql_db.initialize(MySQLDatabase('sndg', user='root', passwd="mito"))

mdb = BioMongoDB("saureus")

BioProjectSubmitter = BioProject.submitters.get_through_model()
BioProjectAssembly = BioProject.assemblies.get_through_model()
# BioProjectAssembly.drop_table()
# BioProjectAssembly.create_table()
esta = 0
noesta = 0

# for x in mdb.db.sequence_collection.find({"ncbi_assembly": {"$exists": 1}}, {"ncbi_assembly": 1}):
#
#     assemblies = list(ExternalAssembly.select().where(ExternalAssembly.assembly_accession == x["ncbi_assembly"]))
#     if assemblies:
#         ea = assemblies[0]
#     else:
#         data = Entrez.read(Entrez.esearch(db="assembly", term=x["ncbi_assembly"] + "[Assembly Accession]"))
#         aid = data["IdList"][0]
#         recurso = Entrez.read(Entrez.esummary(db="assembly", id=aid))
#         data = NCBI().resource_handler["assembly"].attributes(recurso)
#         name = str(NCBI().resource_handler["assembly"].name(recurso))
#         genome = str(data["SpeciesName"])
#         ea = ExternalAssembly(type="assembly", name=name, identifier=aid
#                               , assembly_accession=data['AssemblyAccession'], genome=genome
#                               , assembly_name=data['AssemblyName'])
#         ea.save(force_insert=True)
#
#     links = Entrez.read(Entrez.elink(dbfrom="assembly", id=ea.identifier, linkname="assembly_bioproject"))
#     submitters = []
#     if len(links[0]['LinkSetDb']) > 0:
#         for bioprojlink in links[0]['LinkSetDb']:
#             bioproj_id = bioprojlink["Link"][0]["Id"]
#
#             bioproject = \
#             Entrez.read(Entrez.esummary(db="bioproject", id=bioproj_id))["DocumentSummarySet"]["DocumentSummary"][0]
#             for x in bioproject["Submitter_Organization_List"]:
#                 if list(Submitter.select().where(Submitter.name == x)):
#                     s = Submitter.select().where(Submitter.name == x).get()
#                 else:
#                     s = Submitter(name=x, source="ncbi", rejected=1)
#                     s.save()
#                 submitters.append(s)
#     if not list(ea.submitters):
#         for submitter in submitters:
#             AssemblySubmitters.create(submitter=submitter, resource=ea)


for assembly in ExternalAssembly.select():
    dbassembly = mdb.db.sequence_collection.find_one({"ncbi_assembly": assembly.assembly_accession})
    if dbassembly:
        esta += 1
    else:
        noesta += 1

bioprojects = list(BioProject.select())
for biop in tqdm(bioprojects):

    biop_id = Entrez.read(Entrez.esearch(db="bioproject", term='"' + biop.accession + '"[Accession]'))["IdList"][0]

    biop.identifier = biop_id
    biop.save()

    links = Entrez.read(Entrez.elink(dbfrom="bioproject", id=biop.identifier, linkname="bioproject_assembly_all"))[0][
        'LinkSetDb']
    if links and (not len(list(biop.assemblies))):
        assemblies = []
        for x in links:
            condition = ExternalAssembly.identifier == str(x["Link"][0]["Id"])
            query = list(ExternalAssembly.select().where(condition))
            if query:
                assemblies.append(query[0])
        if assemblies:
            biop.assemblies.add(assemblies)
            biop.save()
    # has_submitter = False
    # for submitter in biop.submitters:
    #     if not submitter.rejected:
    #         has_submitter = True
    #         break
    # if has_submitter:
    #     print biop
