'''
Created on Jun 2, 2014

@author: eze
'''

import json
import logging
import os
import re
import shutil
import subprocess
from tqdm import tqdm
from numpy.lib.scimath import sqrt
from glob import glob
from SNDG import mkdir, execute

_log = logging.getLogger(__name__)

PDB_HEADER_FIELDS = ["Pocket Score",
                     "Drug Score",
                     "Number of V. Vertices",
                     "Mean alpha-sphere radius",
                     "Mean alpha-sphere SA",
                     "Mean B-factor",
                     "Hydrophobicity Score",
                     "Polarity Score",
                     "Volume Score",
                     "Real volume (approximation)",
                     "Charge Score",
                     "Local hydrophobic density Score",
                     "Number of apolar alpha sphere",
                     "Proportion of apolar alpha sphere"]

fpocket_properties_map = {'polar_sasa': 'Polar SASA', 'number_of_alpha_spheres': 'Number of Alpha Spheres',
                          'apolar_alpha_sphere_proportion': 'Apolar alpha sphere proportion',
                          'alpha_sphere_density': 'Alpha sphere density',
                          'charge_score': 'Charge score',
                          'mean_local_hydrophobic_density': 'Mean local hydrophobic density',
                          'total_sasa': 'Total SASA', 'volume': 'Volume',
                          'proportion_of_polar_atoms': 'Proportion of polar atoms',
                          'flexibility': 'Flexibility', 'score': 'Score',
                          'hydrophobicity_score': 'Hydrophobicity score',
                          'apolar_sasa': 'Apolar SASA', 'volume_score': 'Volume score',
                          'cent_of_mass___alpha_sphere_max_dist': 'Cent of mass - Alpha Sphere max dist',
                          'polarity_score': 'Polarity score',
                          'mean_alp_sph_solvent_access': 'Mean alp sph solvent access',
                          'druggability_score': 'Druggability Score',
                          'mean_alpha_sphere_radius': 'Mean alpha sphere radius'}


class Atom(object):
    '''
    Created on Mar 26, 2013
    
    @author: leandro
    '''

    def __init__(self, code, pos, x, y, z, residue=None, bfactor="1.00"):
        '''
        Constructor
        '''
        self.code = code
        self.pos = pos
        self.x = x
        self.y = y
        self.z = z
        self.residue = residue
        self.bfactor = bfactor

    def getInPDBFormat(self, asHetatm=False):
        atom_type = "ATOM  " if not asHetatm else "HETATM"
        atom_pos = "{:5}".format(int(self.pos))
        atom_code = ("   " + self.code)[-3:]
        atom_res_code = ("      " + self.residue.code)[-6:]
        atom_res_chain = ("      " + self.residue.chain)[-2:]
        atom_res_pos = ("      " + str(self.residue.pos))[-4:]
        atom_x = "{:12.3f}".format(self.x)
        atom_y = "{:8.3f}".format(self.y)
        atom_z = "{:8.3f}".format(self.z)
        atom_ocuppancy = "  1.00"
        atom_bfactor = "{:6.2f}".format(float(self.bfactor))
        atom_tail = "           " + self.code[0] + "  "

        return atom_type + atom_pos + atom_code + atom_res_code + atom_res_chain + atom_res_pos + \
               atom_x + atom_y + atom_z + atom_ocuppancy + atom_bfactor + atom_tail

    # @param other: another atom
    # @return: distance between atoms 
    def __sub__(self, other):
        return sqrt(pow(self.x - other.x, 2) + pow(self.y - other.y, 2) + pow(self.z - other.z, 2))

    def __eq__(self, other):
        return self.code == other.code and self.residue == other.residue

    def __repr__(self):
        return "Atom(" + self.code + "," + str(self.pos) + "," + str(self.x) + "," + str(self.y) + "," + str(
            self.z) + "," + str(self.residue) + "," + str(self.bfactor) + ")"

    def __str__(self):
        return "Atom(" + self.code + "," + str(self.pos) + "," + str(self.x) + "," + str(self.y) + "," + str(
            self.z) + "," + str(self.residue) + "," + str(self.bfactor) + ")"


class FpocketOutput():
    def __init__(self, directory):
        self.directory = directory
        self.pdb_id = directory.split("/")[-1].replace("_out", "")
        self.pockets = []
        self.min_pocket_index = 0

    def _info_file_path(self):
        return self.directory + "/{pdb_id}_info.txt".format(pdb_id=self.pdb_id)

    def delete_dir(self):
        shutil.rmtree(self.directory)

    def _parse_info_file(self):
        with open(self._info_file_path()) as f:
            info = f.read()
            info = re.split('\n\s*\n', info)
            for info_pocket in info:
                if info_pocket:
                    lines = info_pocket.split("\n")
                    pocket_num = int(lines[0].split("Pocket")[1].replace(":", "").strip())
                    properties = {prop.split(":")[0].strip().replace(".", ""): prop.split(":")[1].strip() for prop in
                                  lines[1:]}
                    for prop, value in properties.items():
                        try:
                            properties[prop] = float(value)
                        except:
                            pass
                    yield pocket_num, properties

    def parse(self):
        assert os.path.exists(self._info_file_path()), "info file not found in: " + self._info_file_path()

        pockets_index = [int(x.split("_")[0].split("pocket")[1]) for x in os.listdir(self.directory + "/pockets/") if
                         "pocket" in x]
        min_pocket_index = 0
        if pockets_index:
            min_pocket_index = min(pockets_index)

        self.pockets = [p for p in [OutputPocket(pocket_num, self.directory, min_pocket_index, properties) for
                                    pocket_num, properties in self._parse_info_file()] if
                        p.properties["Druggability Score"] > 0.2]

        for pocket in self.pockets:
            pocket.load_atoms()
            pocket.load_alpha()

    @staticmethod
    def load_json(file_path):
        result = FpocketOutput(os.path.dirname(file_path))
        with open(file_path) as h:
            data = json.load(h)
            result.pockets = [p for p in [OutputPocket(pocket_dict["number"], result.directory, 0,
                                                       pocket_dict["properties"], residues=pocket_dict["residues"]
                                                       , atoms=pocket_dict["atoms"]
                                                       ) for pocket_dict in data] if
                              p.properties["Druggability Score"] > 0.2]
        return result

    def save(self, file_path):
        mkdir(os.path.dirname(os.path.abspath(file_path)))
        with open(file_path, "w") as handle:
            json.dump(
                [{"number": p.pocket_num, "residues": p.residues, "as_lines": p.alpha_spheres, "atoms": p.atoms,
                  "properties": p.properties}
                 for p in self.pockets if p.properties["Druggability Score"] > 0.2], handle)

    def __str__(self):
        return "FpocketOutput(%i pockets, max_druggability=%.2f)" % (
            len(self.pockets), max([x.properties["Druggability Score"] for x in self.pockets] or [0]))

    def __repr__(self):
        return self.__str__()


class FPocket(object):
    '''
    Wrapper for fpocket execution
    https://github.com/Discngine/fpocket
    '''

    def __init__(self, pdb_file_path, work_directory=None):
        '''
        Constructor
        '''
        self.pdb_file_path = pdb_file_path
        assert os.path.exists(pdb_file_path), "file %s does not exist" % pdb_file_path
        self.fpocket_binary = "fpocket"
        self.pdb_id = ".".join(os.path.basename(pdb_file_path).split(".")[0:-1])

        self._pdb_file_directory = pdb_file_path.replace(pdb_file_path.split("/")[-1], "")
        if work_directory:
            self.work_directory = work_directory
        else:
            self.work_directory = self._pdb_file_directory
        self.pockets = []

    def _out_directory(self):
        return self.pdb_id + "_out"

    def dest_path(self):
        return self.work_directory + "/" + self._out_directory()

    def hunt_pockets(self):
        abs_path = os.path.abspath(self.pdb_file_path)
        pdb_file = os.path.basename(abs_path)
        pdb_dir = os.path.dirname(abs_path)

        cmd = "docker run -u $(id -u):$(id -g) -w /out -v '{pdb_dir}':/out --rm ezequieljsosa/fpocket {fpocket} -f '{pdb_file}'".format(
            fpocket=self.fpocket_binary, pdb_file=pdb_file, pdb_dir=pdb_dir)
        self._execute(cmd)
        if os.path.abspath(self._pdb_file_directory) != os.path.abspath(self.work_directory):
            if os.path.exists(self.dest_path()):
                shutil.rmtree(self.dest_path(), True)
            work_dir = self._pdb_file_directory + "/" + self._out_directory()
            if os.path.exists(work_dir):
                execute(f'mv "{work_dir}" "{self.dest_path()}"')
        result = FpocketOutput(self.dest_path())
        result.parse()
        return result

    def pdb_out_dir(self):
        return "{dir}/{pdb_id}_out/".format(dir=self.work_directory, pdb_id=self.pdb_id)

    def _execute(self, command):

        _log.debug("Running: " + command)
        try:
            subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True)
        except subprocess.CalledProcessError as e:
            _log.fatal("nCmd error:" + e.output)
            raise
        _log.debug("Command: " + command + "---- Executed correctly")


class OutputPocket(object):
    '''
    classdocs
    '''

    # @param descriptor: list of Residue charges, balls, all the hetatms that defines the pockets
    def __init__(self, pocket_num, directory, min_pocket_index=0, properties=None, residues=None, atoms=None):
        '''
        pocket_num: pocket number in the info file --> Index in the .info and STP file starts in 1
        '''
        self.pocket_num = pocket_num
        self.directory = directory
        self.alpha_spheres = []
        self.atoms = atoms if atoms else []
        self.residues = residues if residues else []
        self.properties = properties if properties else {}

        self.min_pocket_index = min_pocket_index

    def load_atoms(self):
        with open(self._atoms_file_path()) as h:
            for line in h:
                if line[0:6].strip() == 'ATOM':
                    atom_num = line[6:11].strip()
                    self.atoms.append(atom_num)
                    self.residues.append(line[22:26].strip() + line[16].strip())

    def load_alpha(self):
        with open(self._vert_file_path()) as h:
            for line in h:
                if line[0:6].strip() == 'HETATM':
                    # HETATM    1 APOL STP C   1       7.330  72.769  16.106  0.00  0.00          Ve
                    self.alpha_spheres.append(line)
                elif line[0:6].strip() == 'ATOM':
                    # ATOM      1    O STP     1       3.174  33.184  26.211    0.00     4.00
                    line = line.replace("ATOM  ", "HETATM").replace("  O", "POL").replace("   C", "APOL").replace(
                        "  0.00   ", "0.00")
                    self.alpha_spheres.append(line)

    def _atoms_file_path(self):
        return self.directory + "/pockets/pocket{nthpocket}_atm.pdb" \
            .format(nthpocket=self.pocketFileNumber())

    def _out_path(self):
        return glob(self.directory + "/*_out.pdb")[0]

    def _vert_file_path(self):
        return self.directory + "/pockets/pocket{nthpocket}_vert.pqr" \
            .format(nthpocket=self.pocketFileNumber())

    def pocketFileNumber(self):
        '''
        Some versions of fpocket index of the generated files starts the pocket index in 1 and other starts with 0. Ex pocket0.vert.pdb or pocket1.vert.pdb
        '''
        return self.pocket_num - 1 if self.min_pocket_index == 0 else self.pocket_num

    def __str__(self):
        return "OutputPocket(%i pockets, druggability=%.2f)" % (
            self.pocket_num, self.properties["Druggability Score"])

    def __repr__(self):
        return self.__str__()


def fpocket_header_line_handler(header, line):
    if line[0:6] == "HEADER":
        tags = line.split("-")
        if len(tags) > 1:
            prop = line[12:].split(':')
            header[prop[0].strip()] = prop[1].strip()


if __name__ == '__main__':
    from SNDG.PDBs import PDBs

    for pdb in tqdm(PDBs()):
        try:
            pocket_data = "/data/databases/pdb/pockets/" + pdb[1:3] + "/" + pdb + ".json"
            if not os.path.exists(pocket_data):
                fpo = FPocket("/data/databases/pdb/divided/" + pdb[1:3] + "/pdb" + pdb + ".ent")
                res = fpo.hunt_pockets()
                res.save(pocket_data)
        except Exception as e:
            print(e)
