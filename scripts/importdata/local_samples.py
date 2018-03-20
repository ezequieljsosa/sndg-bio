import logging

import xmltodict
from Bio import Entrez
from SNDG import init_log
from SNDG.WebServices.NCBI import ExternalAssembly
from tqdm import tqdm

Entrez.email = "ezejajaja@hotmail.com"

if __name__ == "__main__":
    logger = logging.getLogger('peewee')
    logger.setLevel(logging.INFO)
    init_log()
    search = '("Potentilla"[Organism] OR "Argentina"[Organism] OR argentina[All Fields]) AND "attribute geographic location"[filter]'
    esearch = Entrez.read(Entrez.esearch(db="biosample", term=search, retmax=10000))["IdList"][119:]
    for biosample_id in tqdm(esearch):
        links = Entrez.read(Entrez.elink(dbfrom="biosample", id=biosample_id, linkname="biosample_assembly"))

        if len(links[0]['LinkSetDb']) > 0:
            try:
                biosample = dict(
                    Entrez.read(Entrez.esummary(db="biosample", id=biosample_id))["DocumentSummarySet"]["DocumentSummary"][0])
                acc = biosample["Accession"]
                date = biosample["Date"]
                attributes = dict(xmltodict.parse(biosample["SampleData"]))["BioSample"]["Attributes"]["Attribute"]
                location = [x["#text"] for x in attributes
                            if ("@harmonized_name" in x) and  (x["@harmonized_name"] == "geo_loc_name")][0]
            except KeyError:
                continue
        for assembly_link in links[0]['LinkSetDb']:
            assembly_id = assembly_link["Link"][0]["Id"]
            resource = Entrez.read(Entrez.esummary(db="assembly", id=assembly_id), validate=False)

            data = resource["DocumentSummarySet"]["DocumentSummary"][0]
            name = data["AssemblyName"]
            genome = str(data["SpeciesName"])
            tax = data["Taxid"]
            status = data["AssemblyStatus"]

            exists = ExternalAssembly.select().where(
                ExternalAssembly.assembly_accession == str(data['AssemblyAccession'])).count()

            if not exists:
                ea = ExternalAssembly(type="assembly", name=name, identifier=assembly_id
                                      , assembly_accession=data['AssemblyAccession'], genome=genome
                                      , assembly_name=data['AssemblyName'],
                                      sample_source=acc, sample_date=date, sample_location=location,
                                      assembly_status=status, ncbi_tax=int(tax)
                                      )
                ea.save(force_insert=True)
