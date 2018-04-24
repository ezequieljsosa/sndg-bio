from tqdm import tqdm

from SNDG.Comparative.Pangenome import Pangenome, Strain, sqldb
from SNDG.WebServices.NCBI import NCBI
from Bio import Entrez, SeqIO

from BioSQL import BioSeqDatabase

server = BioSeqDatabase.open_database(driver="MySQLdb", user="root",
                                      passwd="mito", host="localhost", db="bioseqdb")

if __name__ == '__main__':


    from peewee import MySQLDatabase

    mysql_db = MySQLDatabase('bioseqdb', user="root", password="mito")
    sqldb.initialize(mysql_db)
    tables = [Pangenome,Strain]

    # for x in tables:
    #     x.create_table()

    Entrez.email = "Your.Name.Here@example.org"
    query = '"pathogen"[Properties] AND ("Metazoa"[Organism] OR "Viridiplantae"[Organism] OR "Fungi"[Organism] OR "Eukaryota"[Organism] NOT "Metazoa"[Organism] NOT "Fungi"[Organism] NOT "Viridiplantae"[Organism] OR "Bacteria"[Organism] OR txid1224[Orgn] OR "Archaea"[Organism])'
    genomesList = Entrez.read(Entrez.esearch(db="genome", term=query, idtype="acc", retmax=10000))

    genomes = Entrez.read(Entrez.esummary(db="genome", id=",".join(genomesList["IdList"])), validate=False)
    with tqdm(genomes) as pbar1:
        for genome in pbar1:
            # links = Entrez.read(Entrez.elink(dbfrom="genome", id=genome["Id"], linkname="genome_assembly"))
            # links[0]["LinkSetDb"]
            pbar1.set_description(genome["Organism_Name"])
            query = '"%s"[Organism] AND ((latest[filter] OR "latest refseq"[filter]) AND all[filter] NOT anomalous[filter])'
            query = query % genome["Organism_Name"]

            try:
                assembliesList = Entrez.read(Entrez.esearch(db="assembly", term=query, idtype="acc", retmax=100000))
            except:
                with open("errors.log","a") as h:
                    h.write(str(genome["Organism_Name"]) + "\n")
            # {'Number_of_Organelles': '0', 'Number_of_Plasmids': '0', 'Create_Date': '2012/02/29 00:00', 'ProjectID': '0', 'Assembly_Accession': '', 'Number_of_Chromosomes': '1', 'Options': '', 'DefLine': '', u'Item': [], 'Organism_Kingdom': 'Bacteria', 'Assembly_Name': '', 'AssemblyID': '0', 'Organism_Name': 'Mycobacterium paraintracellulare', u'Id': '53912'}
            pan = Pangenome.select().where(Pangenome.name == genome["Organism_Name"])
            if pan.exists():
                pan = pan.get()
            else:
                pan = Pangenome(
                    ncbi_id= int(genome["Id"]),
                    name=genome["Organism_Name"],
                    kingdom=genome["Organism_Kingdom"])
                if "Number_of_Chromosomes" in genome:
                    pan.chromosomes = genome["Number_of_Chromosomes"]
                if "Number_of_Plasmids" in genome:
                    pan.plasmids = genome["Number_of_Plasmids"]
                if "Number_of_Organelles" in genome:
                    pan.organelles = genome["Number_of_Organelles"]
                pan.save()

            with tqdm(assembliesList["IdList"]) as pbar2:
                for assembly_id in pbar2:

                    if Strain.select().where(Strain.ncbi_id == assembly_id).exists():
                        continue

                    try:
                        resource = Entrez.read(Entrez.esummary(db="assembly", id=assembly_id), validate=False)
                    except:
                        with open("errors.log","a") as h:
                            h.write(str(assembly_id) + "\n")
                    # DictElement({u'ChainId': '276825', u'AsmUpdateDate': '2012/07/10 00:00', u'GbUid': '397318', u'PropertyList': ['from-type', 'full-genome-representation', 'has-chromosome', 'latest', 'latest_genbank', 'latest_refseq'], u'SubmissionDate': '2012/02/29 00:00', u'GB_BioProjects': [{u'BioprojectAccn': 'PRJNA82155', u'BioprojectId': '82155'}], u'AssemblyAccession': 'GCF_000276825.1', u'LastMajorReleaseAccession': 'GCF_000276825.1', u'Synonym': {u'RefSeq': 'GCF_000276825.1', u'Genbank': 'GCA_000276825.1', u'Similarity': 'identical'}, u'FtpPath_Regions_rpt': '', u'FtpPath_RefSeq': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/276/825/GCF_000276825.1_ASM27682v1', u'FtpPath_Stats_rpt': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/276/825/GCF_000276825.1_ASM27682v1/GCF_000276825.1_ASM27682v1_assembly_stats.txt', u'AssemblyType': 'haploid', u'AssemblyDescription': '', u'AsmReleaseDate_GenBank': '2012/07/09 00:00', u'ReleaseType': 'Major', u'RS_Projects': [], u'RefSeq_category': 'na', u'LatestAccession': '', u'Coverage': '0', u'AssemblyClass': 'haploid', u'ExclFromRefSeq': [], u'SeqReleaseDate': '2012/02/29 00:00', u'AnomalousList': [], u'AssemblyStatus': 'Complete Genome', u'AsmReleaseDate_RefSeq': '2012/07/10 00:00', u'ReleaseLevel': 'Major', u'Taxid': '1138383', u'PartialGenomeRepresentation': 'false', u'LastUpdateDate': '2012/07/10 00:00', u'RsUid': '398208', u'FromType': 'assembly from type material', u'WGS': '', u'GB_Projects': [], u'BioSampleId': '2603185', u'AssemblyName': 'ASM27682v1', u'EnsemblName': '', u'Organism': 'Mycobacterium paraintracellulare (high GC Gram+)', u'Biosource': {u'Isolate': '', u'InfraspeciesList': [{u'Sub_type': 'strain', u'Sub_value': 'MOTT64'}], u'Sex': ''}, u'SpeciesName': 'Mycobacterium paraintracellulare', u'BioSampleAccn': 'SAMN02603185', u'FtpPath_GenBank': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/276/825/GCA_000276825.1_ASM27682v1', u'ScaffoldN50': '5501090', u'SpeciesTaxid': '1138383', u'UCSCName': '', u'Primary': '398198', u'ContigN50': '5501090', u'Meta': ' <Stats> <Stat category="alt_loci_count" sequence_tag="all">0</Stat> <Stat category="chromosome_count" sequence_tag="all">1</Stat> <Stat category="contig_count" sequence_tag="all">1</Stat> <Stat category="contig_l50" sequence_tag="all">1</Stat> <Stat category="contig_n50" sequence_tag="all">5501090</Stat> <Stat category="non_chromosome_replicon_count" sequence_tag="all">0</Stat> <Stat category="replicon_count" sequence_tag="all">1</Stat> <Stat category="scaffold_count" sequence_tag="all">1</Stat> <Stat category="scaffold_count" sequence_tag="placed">1</Stat> <Stat category="scaffold_count" sequence_tag="unlocalized">0</Stat> <Stat category="scaffold_count" sequence_tag="unplaced">0</Stat> <Stat category="scaffold_l50" sequence_tag="all">1</Stat> <Stat category="scaffold_n50" sequence_tag="all">5501090</Stat> <Stat category="total_length" sequence_tag="all">5501090</Stat> <Stat category="ungapped_length" sequence_tag="all">5501090</Stat> </Stats> <FtpSites>   <FtpPath type="Assembly_rpt">ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/276/825/GCF_000276825.1_ASM27682v1/GCF_000276825.1_ASM27682v1_assembly_report.txt</FtpPath>   <FtpPath type="GenBank">ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/276/825/GCA_000276825.1_ASM27682v1</FtpPath>   <FtpPath type="RefSeq">ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/276/825/GCF_000276825.1_ASM27682v1</FtpPath>   <FtpPath type="Stats_rpt">ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/276/825/GCF_000276825.1_ASM27682v1/GCF_000276825.1_ASM27682v1_assembly_stats.txt</FtpPath> </FtpSites> <assembly-level>5</assembly-level> <assembly-status>Complete Genome</assembly-status> <representative-status>na</representative-status> <submitter-organization>Department of Microbiology and Immunology, Seoul National University</submitter-organization>    ', u'FtpPath_Assembly_rpt': 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/276/825/GCF_000276825.1_ASM27682v1/GCF_000276825.1_ASM27682v1_assembly_report.txt', u'SubmitterOrganization': 'Department of Microbiology and Immunology, Seoul National University', u'RS_BioProjects': [{u'BioprojectAccn': 'PRJNA224116', u'BioprojectId': '224116'}], u'SortOrder': '5C10002768259898'}, attributes={u'uid': u'398208'})
                    data = resource["DocumentSummarySet"]["DocumentSummary"][0]
                    pbar2.set_description(data["AssemblyAccession"])
                    strain = Strain.select().where(Strain.acc == data["AssemblyAccession"])
                    if not strain.exists():
                        strain = Strain(
                            pangenome=pan,
                            acc=data["AssemblyAccession"],
                            tax=data["Taxid"],
                            status=data["AssemblyStatus"],
                            loaded=False,
                            ncbi_id= int( assembly_id)
                        )

                        if "Biosource" in data:
                            if "InfraspeciesList" in data["Biosource"]:
                                strain.name = "|".join( [x["Sub_value"] for x in data["Biosource"]["InfraspeciesList"]] )
                        strain.save()

