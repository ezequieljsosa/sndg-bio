'''
Created on Jun 25, 2015

@author: eze
'''
import logging
import os
import re
import traceback
from multiprocessing.synchronize import Lock
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import CaPPBuilder
from SNDG import Struct
from SNDG import init_log, mkdir, execute
from SNDG.Structure.CompoundTypes import compound_type
from SNDG.Structure.PDBs import PDBs as PDBsIterator
import shutil

init_log()

_log = logging.getLogger("dn_ext")

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

def main(argv=None):  # IGNORE:C0111
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)



    parser = ArgumentParser( formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose", dest="verbose", action="count",
                        help="set verbosity level [default: %(default)s]")

    # parser.add_argument("-dir", "--structs_dir", required = True )
    parser.add_argument("-db", "--database_name", default='pdb')
    parser.add_argument("-host", "--db_host", default='127.0.0.1')

    parser.add_argument( "--procesados", default='/tmp/pdbs_dist_procesados.txt')
    parser.add_argument( "--domains", default='/data/databases/pdb/processed/dns_pdbs.tlb')
    parser.add_argument( "--seqs", default='/data/databases/pdb/processed/pdb_seq_res.fasta')
    parser.add_argument( "--pdbs", default='/data/databases/pdb/')
    parser.add_argument( "--distances", default='/data/databases/pdb/processed/distances.tbl',
                         help="Final output: table with atom distances between residues and ligands. Only for distances less than 'dist' parameter")
    parser.add_argument( "--dist", default=5)
    parser.add_argument( "--pdbs_with_drug", default='/data/databases/pdb/processed/pdbs_with_drug.txt',
                         help="Output: list of PDB codes with an associated ligand")



    args = parser.parse_args()

    if not os.path.exists(args.pdbs):
        sys.stderr.write("%s not found. Specify where is pdbs/divided directory" % (
            parser.pdbs
        ))
        sys.exit(1)
    PDB_PATH = args.pdbs
    CONTACT_DIST = args.dist

    pdbs_with_drug_path = args.pdbs_with_drug
    if not os.path.exists(os.path.dirname(args.pdbs_with_drug)):
        sys.stderr.write("can't %s create %s. Set pdbs_with_drug correctly" % (
            pdbs_with_drug_path
        ))
        sys.exit(1)

    if not os.path.exists(os.path.dirname(args.distances)):
        sys.stderr.write("can't %s create %s. Set distances correctly" % (
            args.distances
        ))
        sys.exit(1)



    pdbs_procesados_path = args.procesados
    print("In %s the processed pdbs are kept, if the file is deleted, the process starts from scratch " % pdbs_procesados_path)
    print("Outputs: '%s' and '%s' " % (pdbs_with_drug_path,args.distances))

    pdbs_procesados = []
    if os.path.exists(pdbs_procesados_path):
        with open(pdbs_procesados_path) as handle:
            pdbs_procesados = [x.strip() for x in handle.readlines()]
        pdbs_procesados = {x: 1 for x in pdbs_procesados}

    pdbs_iterator = PDBsIterator(pdb_dir=args.pdbs)


    def not_processed_iter():
        for pdb, pdb_path in pdbs_iterator:
            if pdb not in pdbs_procesados:
                yield [pdb, pdb_path]

    DNsPDBs = args.domains

    if not os.path.exists(DNsPDBs):
        seqs_from_pdb =args.seqs
        if not os.path.exists(seqs_from_pdb):
            sys.stderr.write("%s does not exists and %s not found. Specify where it is." % (
                DNsPDBs,seqs_from_pdb
            ))
            sys.exit(1)

        sys.stderr.write("%s not found. You can create it with the following command: \n" % DNsPDBs)
        sys.stderr.write(
            "hmmscan --cut_tc --domtblout dns_pdbs.tlb --acc -o pdb_seq_res.hmm Pfam-A.hmm seqs_from_pdb.fasta")
        sys.exit(1)



    drugcompounds = [x for x, y in compound_type.items() if y in ["DRUG", "COFACTOR"]]
    othercompounds = [x for x, y in compound_type.items() if y in ["METAL", "SUGAR", "NUCLEOTIDE", "LIPID"]]
    aminoacidcompounds = [x for x, y in compound_type.items() if y in ["MODIFIED", "RESIDUE"]]

    drugcompounds = othercompounds + drugcompounds

    pdbs_with_drug_path = "/data/databases/pdb/processed/pdbs_with_drug.txt"

    _log.info("proceced pdbs: %i" % len(pdbs_procesados))

    ppb = CaPPBuilder()
    p = PDBParser(PERMISSIVE=1, QUIET=1)

    pdbs_with_drug = []
    if os.path.exists(pdbs_with_drug_path):
        _log.info("pdbs with drugs already loaded")
        with open(pdbs_with_drug_path) as handle:
            for x in handle.readlines():
                pdbs_with_drug.append(x.strip())
    else:
        with open(pdbs_with_drug_path, "a") as handle:
            _log.info("pdbs with drugs will be loaded")
            pdbs = list(pdbs_iterator)
            for pdb, file_path in tqdm(pdbs):
                try:
                    if pdb not in pdbs_with_drug:
                        structure = p.get_structure(pdb, file_path)
                        for res in structure.get_residues():
                            if res.resname in drugcompounds:
                                pdbs_with_drug.append(pdb)
                                handle.write(pdb + "\n")
                                handle.flush()
                                break
                except Exception as ex:
                    print (str(ex))

    # import re
    # dns_table = re.sub(r" +", "\t","\n".join( [str(i) + "\t" + x for i,x in enumerate(open('/data/databases/pdb/processed/dns_pdbs.tlb').readlines()) if not x.startswith("#") ]) )
    if not os.path.exists(DNsPDBs + "2"):
        cols = ["target_name", "accession", "tlen", "query_name", "accession2",
                "qlen", "E-value", "score1", "bias1", "#", "of", "c-Evalue", "i-Evalue",
                "score2", "bias2", "from1", "to1", "from2", "to2", "from3", "to3", "acc"]
        _log.info("correcting hmmer-pdb output")

        regexp = re.compile(" +")
        items = []
        for x in tqdm(open(DNsPDBs).readlines()):
            if not x.startswith("#"):
                line = regexp.split(x)
                items.append(line[0:len(cols)])
                # record = {c: line[i] for i, c in enumerate(cols)}

        df_hmm = pd.DataFrame.from_records(items, columns=cols)
        # df_hmm =  df = pd.read_table('/data/databases/pdb/processed/dns_pdbs.tlb', index_col=None, header=None, delimiter=r"\s+",comment="#",names=cols)
        # df_hmm = df_hmm.dropna()
        df_hmm = df_hmm[["accession", "query_name", "from3", "to3"]]
        df_hmm.to_csv(DNsPDBs + "2")
        df_hmm["pdb"] = map(lambda x: x.split("_")[0].lower().strip(), df_hmm["query_name"])
        df_hmm["chain"] = map(lambda x: x.split("_")[1].upper().strip(), df_hmm["query_name"])
        df_hmm["start_res"] = map(lambda x: x.split("_")[2].upper().strip(), df_hmm["query_name"])
        df_hmm["end_res"] = map(lambda x: x.split("_")[3].upper().strip(), df_hmm["query_name"])
    else:
        df_hmm = pd.read_csv(DNsPDBs + "2")
        df_hmm["pdb"] = map(lambda x: x.split("_")[0].lower().strip(), df_hmm["query_name"])
        df_hmm["chain"] = map(lambda x: x.split("_")[1].upper().strip(), df_hmm["query_name"])
        df_hmm["start_res"] = map(lambda x: x.split("_")[2].upper().strip(), df_hmm["query_name"])
        df_hmm["end_res"] = map(lambda x: x.split("_")[3].upper().strip(), df_hmm["query_name"])
    print (len(df_hmm))

    lock = Lock()


    def centeroid(arr):
        length = len(arr)
        sum_x = np.sum([x.coord[0] for x in arr])
        sum_y = np.sum([x.coord[1] for x in arr])
        sum_z = np.sum([x.coord[2] for x in arr])
        return sum_x / length, sum_y / length, sum_z / length


    def residues_near_drug(drug_centroid, aa_residues):
        residues_near = []
        for r in aa_residues:
            for a in list(r):
                dist = a - Struct(coord=drug_centroid)
                if dist > 20:
                    break
                if dist < 10:
                    residues_near.append(r)
                    break
        return residues_near


    def juan(pdb_raw):
        try:
            pepe(pdb_raw)
        except Exception:
            traceback.print_exc()
        finally:
            with lock:
                pdbs_procesados.append(pdb_raw)
                with open(pdbs_procesados_path, "a") as handle:
                    handle.write(pdb_raw + "\n")


    def pepe(pdb):
        ppb = CaPPBuilder()
        p = PDBParser(PERMISSIVE=1, QUIET=1)
        path_dir = PDB_PATH + "/" + pdb[1:3].lower() + "/"
        path = path_dir + "pdb" + pdb.lower() + ".ent"
        model = list(p.get_structure('X', path))[0]

        for chain_obj in list(model):
            chain = chain_obj.id

            hmm_residues = {}

            pdb_seq = list(model[chain].get_residues())
            if pdb_seq:
                hmm_contacts = {}
                hmm_residues = {}

                hmms = df_hmm[(df_hmm["pdb"] == pdb) & (df_hmm["chain"] == chain) & (
                        df_hmm["start_res"] == str(pdb_seq[0].id[1]))]
                for j, hmm in hmms.iterrows():
                    try:
                        hmm_start = int(hmm["from3"]) - 1
                        hmm_end = int(hmm["to3"]) - 1
                        hmm_chain_name = "_".join(map(str, [hmm["accession"].split(".")[0], hmm["chain"],
                                                            pdb_seq[hmm_start].id[1], pdb_seq[hmm_end].id[1]]))
                        hmm_contacts[hmm_chain_name] = []
                        hmm_residues.update({res.id[1]: hmm_chain_name for res in pdb_seq[hmm_start:hmm_end]})
                    except IndexError:
                        print (pdb, hmm["accession"], hmm["chain"], hmm_start, hmm_end, pdb_seq)

            aa_residues = []
            drug_molecules = []
            for res_obj in chain_obj.get_residues():
                if res_obj.resname in drugcompounds:
                    drug_molecules.append(res_obj)
                elif res_obj.resname in aminoacidcompounds:
                    aa_residues.append(res_obj)

            for res_drug_obj in drug_molecules:
                drug_centroid = centeroid(list(res_drug_obj))
                near_residues = residues_near_drug(drug_centroid, aa_residues)
                for drug_atom in list(res_drug_obj):
                    for near_residue in near_residues:
                        for residue_atom in list(near_residue):
                            distance = (residue_atom - drug_atom)
                            if distance > 20:
                                break
                            if distance < CONTACT_DIST:
                                with open(args.distances, "a") as handle:
                                    hmm_name = hmm_residues[near_residue.id[1]] if near_residue.id[
                                                                                       1] in hmm_residues else "NoDn"
                                    fields = [pdb, chain, hmm_name, near_residue.id[1], near_residue.resname,
                                              residue_atom.serial_number,
                                              res_drug_obj.id[1], res_drug_obj.resname, drug_atom.serial_number, distance]
                                    handle.write("\t".join(map(str, fields)) + "\n")


    _log.info("processing distances file")
    for x in tqdm(set(pdbs_with_drug)):
        if x not in pdbs_procesados:
            juan(x)

    # pool = ThreadPool(1)
    # pool.map(juan, set(pdbs_with_drug) - set(pdbs_procesados))

    print ("Finished!!!")

if __name__ == "__main__":
    sys.exit(main())