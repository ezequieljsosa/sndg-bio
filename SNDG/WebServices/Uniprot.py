import logging
import requests

import Bio.SeqIO as bpio
from io import StringIO
from SNDG.WebServices import download_file
from SNDG import execute

_log = logging.getLogger(__name__)

unip_so_map = {'DNA-binding region': "SO:0001429",
               'active site': "SO:0001104",
               'binding site': "SO:0000409",
               'chain': "SO:0000419",
               'coiled-coil region': "SO:0001080",
               'compositionally biased region': "SO:0001066",
               'cross-link': 'Unip:crosslnk',
               'disulfide bond': "SO:0001088",
               'domain': "SO:0000417",
               'glycosylation site': "Unip:carbohyd",
               'helix': "SO:0001114",
               'initiator methionine': "SO:0000691",
               'lipid moiety-binding region': "Unip:lipid",
               'metal ion-binding site': "SO:0001656",
               'modified residue': "SO:0001089",
               'mutagenesis site': "Unip:mutagen",
               'non-terminal residue': "Unip:non_ter",
               'nucleotide phosphate-binding region': "SO:0001811",
               'propeptide': "SO:0001062",
               'region of interest': "SO:0000001",
               'repeat': "SO:0001068",
               'sequence conflict': "Unip:conflict",
               'sequence variant': "SO:0001060",
               'short sequence motif': "SO:0001067",
               'signal peptide': "SO:0000418",
               'site': "SO:0000408",
               'splice variant': "SO:0001630",
               'strand': "SO:0001111",
               'topological domain': "Unip:topo_dom",
               'transmembrane region': "SO:0001077",
               'turn': "SO:0001128",
               'zinc finger region': "SO:0001971"
               }


class Uniprot(object):
    '''
    classdocs
    '''
    DEFAULT_UNIPROT_URL = 'http://www.uniprot.org/uniprot/'

    def __init__(self):
        '''
        Constructor
        '''
        self.fasta_path = None
        self.output_file = None

        self._query_result = None
        self.queried = False
        self.iterated = False
        self.uniprot_url = Uniprot.DEFAULT_UNIPROT_URL
        self.alias_download_url = "ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/docs/sec_ac.txt"
        self.aliases_begin = 30

    @staticmethod
    def download_fasta(uniprot_id, outdir="./", overwrite=False):
        download_file(Uniprot.DEFAULT_UNIPROT_URL + uniprot_id + ".fasta", f'{outdir}/{uniprot_id}.fasta', overwrite)

    def database_file(self, fasta_path):
        self.fasta_path = fasta_path
        return self

    def download_and_load_seqrecord(self, uniprot_id, format=".xml"):
        res = requests.get(self.uniprot_url + uniprot_id + format, )
        if res.status_code == 200:

                return bpio.read(StringIO(res.text), "uniprot-xml" if "xml" in format else "fasta")

                _log.warn("error parsing: " + uniprot_id)

        return None

    def download_alias(self, dst):
        execute("wget %s -O %s" % (self.alias_download_url, dst), shell=True)

    def load_alias(self, alias_file):
        with open(alias_file, "r") as f:
            lines_demerged = f.readlines()

        self.primaries = {}
        self.dict_accessions = {}
        for l in lines_demerged[self.aliases_begin:]:
            secondary = l.strip().split(" ")[0].strip().lower()
            primary = l.strip().split(" ")[-1].strip().lower()
            if primary not in self.dict_accessions:
                self.dict_accessions[primary] = [primary]
            self.dict_accessions[primary].append(secondary)
            self.primaries[secondary] = primary

    def primary(self, uniprot_id_raw):
        uniprot_id = uniprot_id_raw.lower()
        if uniprot_id in self.primaries:
            return self.primaries[uniprot_id]
        return uniprot_id

    def aliases(self, uniprot_id, include_current=False):
        primary = self.primary(uniprot_id)
        if primary in self.dict_accessions:
            secondaries = self.dict_accessions[primary]
        else:
            secondaries = []
        result = set(secondaries)
        result.add(primary)
        if not include_current:
            if uniprot_id in result:
                result.remove(uniprot_id)
        return result

    @staticmethod
    def download_proteome_from_tax(tax_id, dst_dir, format="fasta"):

        durl = 'http://www.uniprot.org/uniprot/?sort=&desc=&compress=yes&query=taxonomy:{tax}&fil=&format={format}&force=yes'
        download_file(durl.format(tax=tax_id, format=format), dst_dir + "/" + tax_id + "_all.fasta.gz", ovewrite=True)
        execute("gunzip " + dst_dir + "/" + tax_id + "_all.fasta.gz")
        execute("cd-hit -M 0 -c 0.9 -T 0 -i %s -o %s" % (
            dst_dir + "/" + tax_id + "_all.fasta",
            dst_dir + "/" + tax_id + ".fasta"))
        execute("makeblastdb -dbtype prot -in " + dst_dir + "/" + tax_id + ".fasta")

    def blast_para_anotar(self, data_dir, fasta_query, fasta_db):

        execute("makeblastdb -in %s -dbtype prot" % (data_dir + fasta_db))

        blast_result = fasta_db.replace(".fasta", "_blast.xml")

        cmd = "blastp -query %s  -db %s -evalue 0.00001 -outfmt 5  -max_hsps 1 -qcov_hsp_perc 0.8 -num_threads 3 -out %s"
        execute(cmd % (fasta_query, fasta_db, blast_result))
        return blast_result
