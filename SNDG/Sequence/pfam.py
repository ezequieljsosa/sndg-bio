"""
Created on Jun 18, 2014

@author: eze

Sacado de uniprot http://www.ebi.ac.uk/uniprot/TrEMBLstats Release 2016_04 of 13-Apr-2016 (5.  AMINO ACID COMPOSITION  5.1  Composition in percent for the complete database)
"""

import logging
import os
import math
import requests
from requests.exceptions import ConnectionError

from Bio.Data.IUPACData import protein_letters_1to3

abundance = {"A": 8.91, "R": 5.64, "N": 3.98, "D": 5.42, "C": 1.25,
             "E": 6.15, "Q": 3.87, "G": 7.13, "H": 2.21, "I": 5.73,
             "L": 9.85, "K": 5.10, "M": 2.39, "F": 3.95, "P": 4.83, "S": 6.8,
             "T": 5.58, "W": 1.3, "Y": 2.96, "V": 6.8}

_log = logging.getLogger(__name__)


class PFamProfile:

    DEFAULT_TEMPLATE = "/data/databases/xfam/pfam_profiles/%s.hmm"

    @classmethod
    def create_profile(cls, domain_id, hmm_path=None, template=DEFAULT_TEMPLATE):
        domain = domain_id
        if "." in domain_id:
            domain, version = domain_id.split(".")[0:2]  # @UnusedVariable

        hmm_path = hmm_path if hmm_path else template % domain_id
        if not os.path.exists(hmm_path):
            try:
                r = requests.get('http://pfam.xfam.org/family/' + domain + '/hmm')
                if r.status_code == 200:
                    with open(hmm_path, 'w') as f:
                        f.write(r.text)
                else:
                    raise ProfileNotFoundError(domain)
            except ConnectionError:
                raise ProfileNotFoundError(domain)
        return PFamProfile(hmm_path)

    @classmethod
    def split_profiles(cls, pfam_path, out_dir):
        with open(pfam_path) as h:
            txt = h.read()
            profiles = [x.strip() for x in txt.split("//") if x.strip()]
        for p in profiles:
            for x in p.split("\n"):
                if x.startswith("ACC"):
                    acc = x.replace("ACC", "").strip()
                    break
            with open(out_dir + "/" + acc + ".hmm", "w") as h:
                h.write(p + "\n//")

    def __init__(self, profile_path):
        self.profile_path = profile_path
        assert os.path.exists(profile_path)
        self.name = "?"
        self.description = "?"
        self._create_information_dictionary()
        self.abundance = abundance

    def _bitscore_pos(self, pos_aa_dict):
        return sum(
            [math.exp(-x) * math.log(math.exp(-x) / (self.abundance[aa] / 100), 2) for aa, x in pos_aa_dict.items()])

    def bitscores(self):
        return {pos: self._bitscore_pos(aa_dict) for pos, aa_dict in self.info_dict.items()}

    def _entropy(self, pos_aa_dict):
        return - sum([math.exp(-x) * math.log(math.exp(-x), 2) for _, x in pos_aa_dict.items()])

    def important_positions(self, cutoff=1.60):

        return [pos for pos, bitscore in self.bitscores().items() if bitscore > cutoff]

    def _create_information_dictionary(self):
        self.info_dict = {}
        aa_number = 0
        aa_order = sorted(protein_letters_1to3.keys())

        with open(self.profile_path) as handle:
            for line in handle:
                if line.startswith("NAME"):
                    self.name = line.split("NAME")[1].strip()
                if line.startswith("DESC"):
                    self.description = line.split("DESC")[1].strip()
                fields = [x for x in line.split() if x.strip()]  # split('  ')
                try:
                    if line[0:9].strip() == str(aa_number + 1):
                        self.info_dict[aa_number + 1] = dict()
                        actual_field = 1
                        for aa in aa_order:
                            self.info_dict[aa_number + 1][aa] = float(fields[actual_field])
                            actual_field += 1
                        aa_number += 1
                except:
                    logging.warn("errors parsing %s line '%s'" % (self.profile_path, line))

    def items(self, str_alignment):
        aln_items = []
        seq_offset = 0
        profile_offset = -1

        for i, aa in enumerate(str_alignment):

            profile_idx = None
            seq_idx = None

            if aa == ".":  # padding
                seq_offset += 1
                profile_offset += 1
            elif aa == "-":  # deletion
                seq_offset += 1
                profile_idx = i - profile_offset
            elif aa == aa.lower():  # insert
                profile_offset += 1
                seq_idx = i - seq_offset
            else:
                seq_idx = i - seq_offset
                profile_idx = i - profile_offset
                assert aa == aa.upper()  # consensus
            aln_items.append((i, profile_idx, seq_idx, aa))

        return aln_items

    def map_alignment(self, str_alignment, aln_profile_start=0):
        '''http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf.
        In a match column, residues are upper case,
and a '-' character means a deletion relative to the consensus. In an insert column, residues are lower
case, and a '.' is padding. A '-' deletion has a cost: transition probabilities were assessed, penalizing the
transition into and out of a deletion. A '.' pad has no cost per se; instead, the sequence(s) with insertions
are paying transition probabilities into and out of their inserted residue
        '''
        seq_offset = 0
        profile_offset = -1  # pos de profile se numera desde 1
        profile_seq_map = {}
        seq_profile_map = {}
        aln_profile_map = {}

        for i, aa in enumerate(str_alignment):
            if aa == ".":  # padding
                seq_offset += 1
                profile_offset += 1
            elif aa == "-":  # deletion
                seq_offset += 1
                profile_seq_map[i - profile_offset] = None
                aln_profile_map[i] = i - profile_offset + aln_profile_start
            elif aa == aa.lower():  # insert
                profile_offset += 1
            else:
                assert aa == aa.upper()  # consensus
                profile_seq_map[i - profile_offset + aln_profile_start] = i - seq_offset
                seq_profile_map[i - seq_offset] = i - profile_offset + aln_profile_start
                aln_profile_map[i] = i - profile_offset + aln_profile_start

        return profile_seq_map, seq_profile_map, aln_profile_map

    def map_alignment_simple(self, str_alignment, aln_profile_start=0):
        '''
        aln from hmmscan
        '''
        seq_offset = 0
        profile_offset = -1  # pos de profile se numera desde 1
        profile_seq_map = {}
        seq_profile_map = {}
        aln_profile_map = {}

        for i, aa in enumerate(str_alignment):
            if aa == ".":  # padding
                seq_offset += 1
                profile_offset += 1
            elif aa == "-":  # deletion
                seq_offset += 1
                profile_seq_map[i - profile_offset] = None
                aln_profile_map[i] = i - profile_offset + aln_profile_start
            elif aa == aa.lower():  # insert

                seq_profile_map[i - seq_offset] = i - profile_offset + aln_profile_start
                profile_seq_map[i - profile_offset + aln_profile_start] = i - seq_offset
                aln_profile_map[i] = i - profile_offset + aln_profile_start
            else:
                assert aa == aa.upper()  # consensus
                profile_seq_map[i - profile_offset + aln_profile_start] = i - seq_offset
                seq_profile_map[i - seq_offset] = i - profile_offset + aln_profile_start
                aln_profile_map[i] = i - profile_offset + aln_profile_start

        return profile_seq_map, seq_profile_map, aln_profile_map

    def __getitem__(self, key):
        return self.info_dict[key]


class ProfileNotFoundError(Exception):
    pass


if __name__ == '__main__':
    import os
    import argparse

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-p', '--profile_path', default=None)
    parser.add_argument('-n', '--profile_name', default=None)
    parser.add_argument('-t', '--hmm_path_template', default="./%s.hmm")

    args = parser.parse_args()

    assert args.profile_path or args.profile_name, "either profile_name or profile_name must be specified"

    profile = args.profile_name if args.profile_name else args.profile_path.split("/")[-1].split(".")[0]

    pfam_profile = PFamProfile.create_profile(profile, args.profile_path,template=args.hmm_path_template)
    print(pfam_profile.bitscores())
