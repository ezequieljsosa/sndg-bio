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
        self.filter = None

    def make_pdb(self, pdb_path, pdb_id, chain, overwrite=False, struct=None):
        """ Create a new PDB file containing only the specified chains.

        Returns the path to the created file.

        :param pdb_path: full path to the crystal structure
        :param chain: one letter character
        :param overwrite: write over the output file if it exists
        """
        if not self.filter:
            self.filter = SelectChains([chain])

        out_name = "_".join([pdb_id, chain + ".pdb"])
        out_path = os.path.join(self.out_dir, out_name)
        _log.debug("OUT PATH:" + out_path)

        # Skip PDB generation if the file already exists
        if (not overwrite) and (os.path.isfile(out_path)):
            _log.debug("Chain %s of '%s' already extracted to '%s'." %
                       (chain, pdb_id, out_name))
            return out_path

        _log.debug("Extracting chain %s from %s..." % (
            chain, pdb_id))

        # Get structure, write new file with only given chains
        if struct is None:
            struct = self.parser.get_structure(pdb_id, pdb_path)
        self.writer.set_structure(struct)
        self.writer.save(out_path, select=self.filter)

        return out_path


if __name__ == "__main__":
    from SNDG import init_log

    init_log()

    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--inputpdb", required=True)
    parser.add_argument("-c", "--chain", required=True)
    parser.add_argument("-o", "--outdir", default="./")

    args = parser.parse_args()

    pdb = ".".join( args.inputpdb.split("/")[-1].split(".")[:-1])

    splitter = ChainSplitter(args.outdir)

    splitter.make_pdb(args.inputpdb,pdb,chain=args.chain,overwrite=True)
