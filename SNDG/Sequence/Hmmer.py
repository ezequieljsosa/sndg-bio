'''
Created on Mar 17, 2014

@author: eze
'''

import logging
import os
import subprocess

import Bio.SeqIO as bpio
from Bio.SearchIO.HmmerIO.hmmer3_text import Hmmer3TextParser
from Bio.SeqRecord import SeqRecord
from SNDG import Struct
from SNDG.Sequence.so import SO_TERMS

_log = logging.getLogger(__name__)

DEFAULT_PARAMS = {"--acc": None, "--cut_tc": None, "--notextw": None, "--cpu": "4"}
DEFAULT_PATHS = {"hmmscan": "hmmscan", "hmmpress": "hmmpress"}
DEFAULT_OUTPUT_PATH = "hmmresult.hmm"

HMM_BINARIES_EXTENTIONS = [".h3f", ".h3i", ".h3m", ".h3p"]


class Hmmer(object):
    '''

    
    '''

    def __init__(self, query, database="/data/databases/xfam/Pfam-A.hmm", output_file=None, params=DEFAULT_PARAMS,
                 binaries=DEFAULT_PATHS):
        '''

        '''
        self.database = database
        self.seq_query = query
        self.params = params
        self.binary_paths = binaries

        self.output_file = output_file

        #             self.output_file = self.pipeline.tmpfile(".hmm")
        self.so_type = SO_TERMS["polypeptide_domain"]
        self.cmd = None
        # self.query_iterator = None
        self.query_iterator = self._parse_output

    def _check_input_file(self):
        if not os.path.exists(self.database):
            msg = "file %s doest not exists" % self.database
            _log.error(msg)
            raise IOError(msg)
        if not os.path.exists(self.seq_query):
            msg = "file %s doest not exists" % self.seq_query
            _log.error(msg)
            raise IOError(msg)

    def _format_database(self):
        if not all(map(os.path.exists, self._database_binary_files())):
            self.hmmpress()

    def _database_binary_files(self):
        return [self.database + x for x in HMM_BINARIES_EXTENTIONS]

    def create_domains_fasta(self, domains_fasta_path):

        if (not os.path.exists(domains_fasta_path)) or (os.path.getsize(domains_fasta_path) == 0):
            sequences = {x.id: x.seq for x in bpio.parse(self.seq_query, "fasta")}
            with open(domains_fasta_path, "w") as handle:
                for query in self._parse_output():
                    for hit in query:
                        for hsp in hit:
                            dn_id = "_".join([query.id, hit.id, str(hsp.query_start), str(hsp.query_end)])
                            record = SeqRecord(seq=sequences[query.id][hsp.query_start:hsp.query_end], id=dn_id)
                            bpio.write(record, handle, "fasta")
            _log.debug("dominios fueron generados en %s" % domains_fasta_path)

    def hmmpress(self):
        '''
        Run example
            hmmpress pfamA.hmm
        '''
        command = "{binary} '{database}' " \
            .format(binary=self.binary_paths["hmmpress"], database=self.database)
        self._execute(command)

    def query(self):
        "hmmscan [-options] <hmmdb> <seqfile>"
        self._format_database()
        self._check_input_file()

        if (not os.path.exists(self.output_file)) or (os.path.getsize(self.output_file) == 0):
            self._hmmscan()

        decorated = self.query_iterator()
        meta = self.query_iterator()._meta

        return Struct(_meta=meta, __iter__=lambda: decorated)

    def _hmmscan(self):
        '''
        Run example
            hmmscan --acc --cut_tc --noali  pfamA.hmm arch_id.fasta
        '''
        self.cmd = "{binary} {params} -o '{output_file}' '{database}'  '{aa_input}'" \
            .format(binary=self.binary_paths["hmmscan"], database=self.database,
                    params=self._params2str(), output_file=self.output_file, aa_input=self.seq_query)
        self._execute(self.cmd)
        # TODO: hacer esto mismo con python
        command = "grep -v '{str_remove}' '{output_file}' > '{output_file}2';mv '{output_file}2' '{output_file}'".format(
            str_remove='\[ok\]', output_file=self.output_file)
        self._execute(command)

        return self._parse_output()

    def _complete_data(self, query_result):
        '''
        Standarizes the fields that the objects should have:
        Query - accession = SO:0000417 
        '''
        # for query_result in query_map.values():
        for hit in query_result:
            hit.accession = SO_TERMS["polypeptide_domain"]
        return query_result

    def _execute(self, command):

        _log.debug("Running: " + command)
        try:
            subprocess.call(command, stderr=subprocess.STDOUT, shell=True, executable='/bin/bash')
        except subprocess.CalledProcessError as e:
            _log.fatal("nCmd error:" + e.output)
            raise
        _log.debug("Command: " + command + "---- Executed correctly")

    def _parse_output(self):
        # SearchIO.parse(self.output_file, 'hmmer3-text')
        output_file = open(self.output_file)
        return Hmmer3TextParser(output_file)

    def clean_temp(self):
        try:
            os.remove(self.output_file)
        except Exception as ex:
            _log.warn(ex)

    def clean_db(self):
        try:
            for file2remove in self._database_binary_files():
                os.remove(file2remove)
        except Exception as ex:
            _log.warn(ex)

    def _params2str(self):
        str_params = ""
        for key, val in self.params.items():
            if val:
                str_params += key + "=" + val
            else:
                str_params += key
            str_params += " "
        return str_params
