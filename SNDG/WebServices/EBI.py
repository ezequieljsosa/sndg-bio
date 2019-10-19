"""

"""
from tqdm import tqdm
from SNDG.WebServices import download_file
import requests
import os
import logging
import subprocess as sp

_log = logging.getLogger(__name__)


class EBI:

    # sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp&download

    @staticmethod
    def download_ena_project(project_id, dst_dir):
        dst_dir = os.path.abspath(dst_dir)
        url_template = "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=" + project_id + "&result=read_run&fields=sample_accession,experiment_accession,run_accession,fastq_ftp,fastq_md5&download=txt"
        r = requests.get(url_template)
        if r.status_code == 200:
            lines = r.text.split("\n")
            with tqdm(lines) as pbar:
                for l in pbar:

                    if len(l.strip().split("\t")) > 3:
                        sample_accession, experiment_accession, run_accession, fastq_ftp,fastq_md5 = l.split("\t")
                        if len(fastq_ftp.split(";")) == 2:
                            basefilename = dst_dir + "/" + "_".join([sample_accession, experiment_accession, run_accession])

                            if (not os.path.exists(basefilename + "_1.fastq.gz")) and (not os.path.exists(basefilename + "_1.fastq")):
                                pbar.set_description(fastq_ftp.split(";")[0])
                                try:
                                    download_file(fastq_ftp.split(";")[0],
                                                  basefilename + "_1.fastq.gz")
                                except:
                                    _log.warn("error downloading: " + basefilename + "_1.fastq.gz")
                                    try:
                                        os.rmdir(basefilename + "_1.fastq.gz")
                                    except:
                                        pass
                            f1md5 = sp.check_output("md5sum %s_1.fastq.gz" % basefilename,shell=True).split()[0].strip()
                            if fastq_md5.split(";")[0] != f1md5:
                                print("%s error md5 sum" % basefilename)
                                try:
                                    os.rmdir(basefilename + "_1.fastq.gz")
                                except:
                                    pass

                            if (not os.path.exists(basefilename + "_2.fastq.gz")) and (not os.path.exists(basefilename + "_2.fastq")):
                                pbar.set_description(fastq_ftp.split(";")[1])
                                try:
                                    download_file(fastq_ftp.split(";")[1],
                                                  basefilename + "_2.fastq.gz")
                                except:
                                    _log.warn("error downloading: " + basefilename + "_2.fastq.gz")
                                    try:
                                        os.rmdir(basefilename + "_2.fastq.gz")
                                    except:
                                        pass
                            f1md5 = sp.check_output("md5sum %s_2.fastq.gz" % basefilename,shell=True).split()[0].strip()
                            if fastq_md5.split(";")[1] != f1md5:
                                print("%s error md5 sum" % basefilename)
                                try:
                                    os.rmdir(basefilename + "_2.fastq.gz")
                                except:
                                    pass

        else:
            raise Exception("request error %i" % r.status_code)


if __name__ == '__main__':
    EBI.download_ena_project("PRJEB7669","./")
