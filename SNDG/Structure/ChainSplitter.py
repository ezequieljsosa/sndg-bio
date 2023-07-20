'''
Created on Aug 1, 2016

https://stackoverflow.com/questions/11685716/how-to-extract-chains-from-a-pdb-file

@author: David Cain
'''

import os
import logging
from argparse import ArgumentParser, RawDescriptionHelpFormatter

from Bio import PDB
from Bio.Data.IUPACData import protein_letters_1to3

_log = logging.getLogger(__name__)


class SelectChains(PDB.Select):
    """ Only accept the specified chains when saving. """

    def __init__(self, chain_letters):
        self.chain_letters = chain_letters
        self.valid_resnames = [x.upper() for x in protein_letters_1to3.values()]

    def accept_chain(self, chain):
        return (chain.get_id() in self.chain_letters)

    def accept_residue(self, residue):
        return ((not bool(residue.id[0].strip()))
                and (residue.resname in self.valid_resnames)
                )


class SelectResidues(PDB.Select):
    """ Only accept the specified chains when saving. """

    def __init__(self, chain_letters, resids):
        self.chain_letters = chain_letters
        self.resids = resids
        self.valid_resnames = [x.upper() for x in protein_letters_1to3.values()]

    def accept_chain(self, chain):
        return (chain.get_id() in self.chain_letters)

    def accept_residue(self, residue):
        return ((not bool(residue.id[0].strip()))
                and (residue.resname in self.valid_resnames)
                and residue.id[1] in self.resids
                )


class ChainSplitter:
    def __init__(self, out_dir=None):
        """ Create parsing and writing objects, specify output directory. """
        self.parser = PDB.PDBParser(PERMISSIVE=True, QUIET=True)
        self.writer = PDB.PDBIO()
        if out_dir is None:
            out_dir = os.path.join(os.getcwd(), "chain_PDBs")
        self.out_dir = out_dir

    def extract_chain(self, pdb_id, chain, out_path, pdbfilter, overwrite=False, struct=None):

        _log.debug("Extracting chain %s from %s..." % (
            chain, pdb_id))

        # Get structure, write new file with only given chains

        self.writer.set_structure(struct)
        self.writer.save(out_path, select=pdbfilter)

        return out_path

    def make_pdb(self, pdb_path, pdb_id, chains=None, overwrite=False, struct=None):
        """ Create a new PDB file containing only the specified chains.

        Returns the path to the created file.

        :param pdb_path: full path to the crystal structure
        :param chain: one letter character
        :param overwrite: write over the output file if it exists
        """

        if struct is None:
            struct = self.parser.get_structure(pdb_id, pdb_path)
        if chains:
            chains = chains.split(",")
        else:
            chains = [c.id for c in list(struct[0])]

        for c in chains:
            out_name = "_".join([pdb_id, c + ".pdb"])
            out_path = os.path.join(self.out_dir, out_name)


            # Skip PDB generation if the file already exists
            if (not overwrite) and (os.path.isfile(out_path)):
                _log.debug("Chain %s of '%s' already extracted to '%s'." %
                           (c, pdb_id, out_name))
                continue
            _log.debug("OUT PATH:" + out_path)
            pdbfilter = SelectChains(c)
            self.extract_chain(pdb_id, c, out_path, pdbfilter, overwrite, struct)


if __name__ == "__main__":
    from SNDG import init_log, mkdir

    init_log()

    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("inputpdb")
    parser.add_argument("-c", "--chains", required=False, default="",
                        help="coma separated list of chains")
    parser.add_argument("-o", "--outdir", default="./")
    parser.add_argument("--overwrite", action="store_true")


    args = parser.parse_args()

    mkdir(args.outdir)
    assert os.path.exists(args.outdir), f'"{args.outdir}" could not be created'
    pdb = ".".join(args.inputpdb.split("/")[-1].split(".")[:-1])

    splitter = ChainSplitter(args.outdir)

    splitter.make_pdb(args.inputpdb, pdb, chains=args.chains, overwrite=args.overwrite)
