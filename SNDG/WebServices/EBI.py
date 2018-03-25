"""

"""
from tqdm import tqdm
from SNDG.WebServices import download_file
import requests


class EBI:

    # sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp&download

    @staticmethod
    def download_ena_project(project_id):
        url_template = "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=" + project_id + "&result=read_run&fields=sample_accession,experiment_accession,run_accession,fastq_ftp&download=txt"
        r = requests.get(url_template)
        if r.status_code == 200:
            lines = r.text.split("\n")
            for l in tqdm(lines):
                if len(l.strip().split("\t")) > 3:
                    sample_accession, experiment_accession, run_accession, fastq_ftp = l.split("\t")
                    if len(fastq_ftp.split(";")) == 2:
                        basefilename = "_".join([sample_accession, experiment_accession, run_accession])
                        download_file(fastq_ftp.split(";")[0],
                                      basefilename + "_1.fastq.gz")
                        download_file(fastq_ftp.split(";")[1],
                                      basefilename + "_2.fastq.gz")

        else:
            raise Exception("request error %i" % r.status_code)


if __name__ == '__main__':
    EBI.download_ena_project("PRJEB7669")
