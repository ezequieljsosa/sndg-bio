"""

"""
from SNDG import mkdir,init_log
from SNDG.WebServices import download_file


class EBI:

    @staticmethod
    def download_ena_project(project_id):
        url_templete = "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=" + project_id + "&result=read_run&fields=sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp&download=txt"
