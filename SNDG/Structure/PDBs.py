'''
Created on Jul 4, 2014

@author: eze
'''
import logging
import os
from collections import defaultdict

import Bio.SeqIO as bpio
import pandas as pd
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import PDBList
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import is_aa, CaPPBuilder
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import seq1
from tqdm import tqdm

from SNDG import mkdir, execute
from SNDG.WebServices import download_file
import fileinput

from SNDG.Structure.ChainSplitter import ChainSplitter

_log = logging.getLogger(__name__)

SIFTS_GO_TERMS = {"loop": "SO:0100012", "strand": "SO:0001111", "helix": "SO:0001114"}
HOMOLOGOUS = ["InterPro", "UniProt", "SCOP", "CATH", "NCBI"]


# ftp://ftp.ebi.ac.uk/pub/databases/pdb/derived_data/index/entries.idx

class PDBFile:

    def __init__(self, code, in_pdb):
        self.code = code
        self.struct = PDBParser(PERMISSIVE=1, QUIET=1).get_structure(self.code, in_pdb)

    def residues_map(self, selected_chain=None, standard_aa=True):
        rmap = {}
        for chain in self.struct.get_chains():
            if (not selected_chain) or (selected_chain == chain.id):
                residues = [x for x in chain.get_residues() if is_aa(x, standard=standard_aa)]
                rmap[chain.id] = {i: x.id for i, x in enumerate(residues)}
        return rmap

    def seq(self, selected_chain=None, standard_aa=True):
        records = []
        for chain in self.struct.get_chains():
            if (not selected_chain) or (selected_chain == chain.id):
                residues = [x for x in chain.get_residues() if is_aa(x, standard=standard_aa)]
                if residues:
                    seq = "".join([seq1(x.resname) for x in residues])
                    start = str(residues[0].id[1])
                    end = str(residues[-1].id[1])
                    record = SeqRecord(id="_".join([self.code, chain.id, start, end]), description="", seq=Seq(seq))
                    records.append(record)
        return records


class PDBs(object):
    '''

    '''

    def __iter__(self):
        for index_dir in os.listdir(self.pdbs_dir):
            if len(index_dir) == 2:
                for x in os.listdir(self.pdbs_dir + "/" + index_dir):
                    if x.endswith(self.pdb_extention):
                        yield x[3:7], self.pdbs_dir + "/" + index_dir + "/" + x

    def __init__(self, pdb_dir='/data/databases/pdb/'):
        '''

        '''
        self.pdb_dir = pdb_dir + ("/" if pdb_dir[-1] != "/" else "")
        self.pdbs_dir = self.pdb_dir + 'divided/'
        self.pockets_dir = self.pdb_dir + 'pockets/'

        self.pdb_seq_res_path = self.pdb_dir + "/pdb_seqres.txt"
        self.url_pdb_seq_res = "ftp://ftp.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt"
        self.url_pdb_entries = 'ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/entries.idx'

        self.url_pdb_files = 'ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/'  # 00/pdb100d.ent.gz

        self.pdb_extention = ".ent"
        self.pdb_download_extention = ".ent.gz"

        self.uncompress_file = True
        self.delete_compressed = True
        self.entries_path = pdb_dir + '/entries.idx'
        self._entries_df = None

    def entries_df(self):
        if not hasattr(self._entries_df, "apply"):
            entries_columns = ["IDCODE", "HEADER", "ACCESSIONDATE", "COMPOUND", "SOURCE", "AUTHORS", "RESOLUTION",
                               "EXPERIMENT"]
            data = []
            with open(self.entries_path) as h:
                for i, line in enumerate(h):
                    if i < 3:
                        continue
                    vec = line.split("\t")
                    data.append({c: vec[idx] for idx, c in enumerate(entries_columns)})

            # self._entries_df = pd.read_table(self.entries_path, skiprows=[0, 1, 2], sep='\t', names=entries_columns)
            self._entries_df = pd.DataFrame(data)
        return self._entries_df

    def pdb_path_base(self, pdb):
        return self.pdbs_dir + "/" + pdb[1:3]

    def pdb_path(self, pdb):
        return self.pdbs_dir + "/" + pdb[1:3] + "/pdb" + pdb + self.pdb_extention

    def pdb_path_gzipped(self, pdb):
        return self.pdb_path(pdb) + ".gz"

    def pdb_pockets_path(self, pdb):
        return self.pockets_dir + "/" + pdb[1:3] + "/pdb" + pdb + self.pdb_extention + ".json"

    def download_pdb_seq_ses(self):
        download_file(self.url_pdb_seq_res, self.pdb_seq_res_path, ovewrite=True)

    def download_pdb_entries(self):
        download_file(self.url_pdb_entries, self.entries_path, ovewrite=True)

    def update_pdb_dir(self):
        assert os.path.exists(self.pdb_dir), "the target directory does not exists %s" % self.pdb_dir
        assert os.path.exists(self.entries_path), "the entries path does not exists %s" % self.entries_path

        pdbs = [x.lower() for x in self.entries_df().IDCODE]
        pbar = tqdm(pdbs)
        for pdb in pbar:
            try:
                pbar.set_description(pdb)
                self.update_pdb(pdb)

            except Exception as ex:
                _log.warn(str(ex))

    def update_pdb(self, pdb):
        pdb = pdb.lower()
        mkdir(self.pdbs_dir + pdb[1:3])
        if not os.path.exists(self.pdb_path(pdb)) or (os.path.getsize(self.pdb_path(pdb)) < 100):
            if os.path.exists(self.pdb_path_gzipped(pdb)) and (os.path.getsize(self.pdb_path_gzipped(pdb)) > 100):
                execute("gunzip " + self.pdb_path_gzipped(pdb))
                if os.path.exists(self.pdb_path_gzipped(pdb)) and not os.path.exists(self.pdb_path(pdb)):
                    os.remove(self.pdb_path_gzipped(pdb))
            elif not os.path.exists(self.pdb_path(pdb)):
                download_file(self.url_pdb_files + pdb[1:3] + "/pdb" + pdb + self.pdb_download_extention,
                              self.pdbs_dir + pdb[1:3] + "/pdb" + pdb + self.pdb_download_extention, ovewrite=True)
                execute("gunzip " + self.pdb_path_gzipped(pdb))
        return self.pdb_path(pdb)

    def pdbs_seq_for_modelling(self, out_fasta=None,
                               pdbsIter=None, reuse_previours=None):
        if pdbsIter == None:
            pdbsIter = PDBs(self.pdb_dir)
        if not out_fasta:
            out_fasta = self.pdb_dir + "processed/seqs_from_pdb.fasta"

        pdb_codes = {x.lower(): 1 for x in self.entries_df().IDCODE}

        reuse = defaultdict(lambda: [])
        if reuse_previours:
            for x in bpio.parse(reuse_previours, "fasta"):
                pdb = x.id.split("_")[0]
                reuse[pdb].append(x)
        reuse = dict(reuse)

        pdblist = list(pdbsIter)
        with open(out_fasta, "w") as out_fasta_handle:
            for (pdb, pdb_file_path) in tqdm(pdblist):
                if pdb in pdb_codes:
                    if pdb in reuse:
                        bpio.write(reuse[pdb], out_fasta_handle, "fasta")
                    else:
                        self.seq_from_pdb(out_fasta_handle, pdb, pdb_file_path)

    def records_from_pdb(self, pdb, pdb_file_path, standard_aa=True, selected_chain=None):
        records = []
        struct = PDBParser(PERMISSIVE=1, QUIET=1).get_structure(pdb, pdb_file_path)
        for chain in struct.get_chains():
            if (not selected_chain) or (selected_chain == chain.id):
                residues = [x for x in chain.get_residues() if is_aa(x, standard=standard_aa)]
                if residues:
                    seq = "".join([seq1(x.resname) for x in residues])
                    start = str(residues[0].id[1])
                    end = str(residues[-1].id[1])
                    record = SeqRecord(id="_".join([pdb, chain.id, start, end]), description="", seq=Seq(seq))
                    records.append(record)
        return records

    def seq_from_pdb(self, out_fasta_handle, pdb, pdb_file_path, standard_aa=True, selected_chain=None):

        records = self.records_from_pdb(pdb, pdb_file_path, standard_aa, selected_chain)
        for record in records:
            bpio.write(record, out_fasta_handle, "fasta")

    def extract_chains(self, pdb_chain_list, output_dir):
        for pdb, chain in pdb_chain_list:
            ChainSplitter(output_dir).make_pdb(self.pdb_path(pdb), pdb, chain, overwrite=True)

    def pdb_path(self, pdb):
        return self.pdbs_dir + pdb[1:3] + "/pdb" + pdb + ".ent"

    @staticmethod
    def sequence_from_residues(residues):
        return "".join([protein_letters_3to1[res.get_resname()[0] + res.get_resname()[1:3].lower()]
                        for res in residues])


if __name__ == '__main__':
    import sys
    from SNDG import init_log

    init_log()

    import argparse
    from SNDG.Structure.PDBs import PDBs

    parser = argparse.ArgumentParser(description='PDB utils')

    subparsers = parser.add_subparsers(help='commands', description='valid subcommands', dest='command')

    update_pdb = subparsers.add_parser('update', help='List contents')
    update_pdb.add_argument('-i', '--pdbs_dir', help="pdbs_directory", default="/data/databases/pdb/")

    args = parser.parse_args()
    if args.command == "update":
        # remzemeber to configure ftp
        pdbs = PDBs(pdb_dir=args.pdbs_dir)
        pdbs.download_pdb_entries()
        pdbs.update_pdb_dir()
        pdbs.download_pdb_seq_ses()
        sys.exit(0)

    # os.environ["ftp_proxy"] = "http://proxy.fcen.uba.ar:8080"
    # pdbs.download_pdb_seq_ses()

    # from SNDG.Structure.PDBs import PDBs
    # pdbs = PDBs(pdb_dir="/data/databases/pdb/")
    # pdbs.pdbs_seq_for_modelling("/data/databases/pdb/processed/seqs_from_pdb.fasta")
    # pepe = pdbs.entries_df()
    # print pepe
