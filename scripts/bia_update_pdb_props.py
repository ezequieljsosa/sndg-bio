#!/usr/bin/python
# encoding: utf-8
'''
scripts.load_bia_proteome -- shortdesc

scripts.load_bia_proteome is a description

It defines classes_and_methods

@author:     Ezequiel Sosa
@copyright:  2015 BIA. All rights reserved.
@license:    license
@contact:    user_email

'''
import logging
import os
import sys
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import pandas as pd
from Bio.PDB.PDBExceptions import PDBConstructionException
from tqdm import tqdm
import Bio.SearchIO as bpsio

# from BiaStructure.IO.pdb import PDBs
from Bio.PDB.PDBParser import PDBParser
from mongoengine.errors import DoesNotExist

from SNDG import init_log
from SNDG.BioMongo.Model import BioProperties
from SNDG.BioMongo.Model import Cluster
from SNDG.BioMongo.Model.Structure import ExperimentalStructure, ResidueSet
from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
from SNDG.Structure.CompoundTypes import get_compound_type
from SNDG.Structure.FPocket import FPocket
from SNDG.Structure.PDBs import PDBs


# from Bia.Programs.Cluster.CDHit import CDHit

init_log()
_log = logging.getLogger(__name__)

from SNDG.Structure.CompoundTypes import compound_type


def update_clusters():
    for cluster_name, seqs in CDHit().clustered_seq_iterator("/data/databases/pdb/processed/seqs_from_pdb95.fasta"):
        _log.debug(cluster_name)

        cristals = []
        cluster = Cluster(name=cluster_name, type="PDB_Segments_95")
        for seq in seqs:
            seq_id, seq_start, seq_end, clust_start, clust_end = seq
            pdb, chain, start, end = seq_id.split("_")
            cristals.append(pdb)
            cluster.parts.append(BioProperties(pdb=pdb, chain=chain, start=start, end=end,
                                               seq_start=seq_start, seq_end=seq_end, clust_start=clust_start,
                                               clust_end=clust_end))
        for pdb in set(cristals):
            try:
                cristal_doc = ExperimentalStructure.objects(name=pdb).no_cache().get()
                cristal_doc.clusters = [x for x in cristal_doc.clusters if x.type != "PDB_Segments_95"]
                if not cristal_doc.cluster(cluster_name):
                    cristal_doc.clusters.append(cluster)
                    cristal_doc.save()
            except DoesNotExist as ex:
                print(str(ex))


def update_binding_residues(distances_tbl):
    df_binding_dist = pd.read_table(distances_tbl, sep="\t",
                                    names=[
                                        "pdb", "chain", "hmm_name", "prot_res", "resname",
                                        "res_atom_id",
                                        "comp_res_id", "comp_resname", "comp_atom_id", "distance"
                                    ])
    df_binding_dist["comptype"] = map(lambda x: compound_type[x], df_binding_dist.comp_resname)

    groups = df_binding_dist.groupby("pdb")
    total = len(groups)
    _log.debug("procesing binding")
    for pdb, df_binding_dist_pdb in tqdm(groups, total=total):

        try:
            for r_comp_type in ['LIPID', 'METAL', 'NUCLEOTIDE', 'SUGAR', "DRUG", "COFACTOR"]:
                if r_comp_type != "COFACTOR":
                    comp_type = r_comp_type
                else:
                    comp_type = "DRUG"

                df_comp_dist_pdb_near = df_binding_dist_pdb[
                    (df_binding_dist_pdb.distance <= 3) & (comp_type == df_binding_dist_pdb.comptype)]
                if len(df_comp_dist_pdb_near):
                    strdoc = ExperimentalStructure.objects(name=pdb).no_cache().get()
                    rs_name = comp_type.lower() + "_binding"
                    if not strdoc.has_residue_set(rs_name):
                        residue_list = list(set([row.chain + "_" + str(row.prot_res) for i, row in
                                                 df_comp_dist_pdb_near.iterrows()]))  # @UnusedVariable
                        rs = ResidueSet(name=rs_name, residues=residue_list, type="binding")
                        strdoc.residue_sets.append(rs)
                        strdoc.save()

        except DoesNotExist:
            _log.warn("%s does not exists" % pdb)


def free_cys_tyr(pdb_utils):
    parser = PDBParser(PERMISSIVE=1, QUIET=1)
    _log.debug("procesing free cys/tyr")
    total = ExperimentalStructure.objects(residue_sets__name__ne = "free_tyr").count()
    for strdoc in tqdm(ExperimentalStructure.objects(residue_sets__name__ne = "free_tyr").no_cache().timeout(False), total=total):

        if not (strdoc.residue_set("free_cys") or strdoc.residue_set("free_tyr")):
            if not os.path.exists(pdb_utils.pdb_path(strdoc.name)):
                pdb_utils.update_pdb(strdoc.name)
            if not os.path.exists(pdb_utils.pdb_path(strdoc.name)):
                continue
            try:
                bp_pdb = list(parser.get_structure(strdoc.name, pdb_utils.pdb_path(strdoc.name)  ))[0]
            except PDBConstructionException:
                continue
            except TypeError:
                continue

            free = {"CYS": [], "TYR": []}
            codes = {"CYS": "SG", "TYR": "OH"}
            for x in bp_pdb.get_residues():
                if x.resname in codes:
                    neighbor_atoms = set(list(bp_pdb.get_atoms())) - set(list(x))
                    if (codes[x.resname] in x) and (
                            not any(map(lambda atom: (x[codes[x.resname]] - atom) <= 3, neighbor_atoms))):
                        free[x.resname].append(x.parent.id + "_" + str(x.id[1]))
            if free["CYS"]:
                rs = ResidueSet(name="free_cys", residues=free["CYS"])
                strdoc.residue_sets.append(rs)
            if free["TYR"]:
                rs = ResidueSet(name="free_tyr", residues=free["TYR"])
                strdoc.residue_sets.append(rs)
            if free["CYS"] or free["TYR"]:
                strdoc.save()


def update_csa(csa_txt):
    df_csa = pd.read_csv(csa_txt)

    for pdb in tqdm(set(df_csa["PDB ID"])):
        try:
            strdoc = ExperimentalStructure.objects(name=pdb).no_cache().get()
        except ExperimentalStructure.DoesNotExist:
            _log.warn(pdb + " csa pdb does not exists...")
            continue

        pdb = strdoc.name
        pdb_csa = df_csa[df_csa["PDB ID"] == pdb]
        if len(pdb_csa) > 0:
            if not strdoc.has_residue_set("csa"):
                csas = [str(row["CHAIN ID"]) + "_" + str(row["RESIDUE NUMBER"]) for i, row in
                        pdb_csa.iterrows()]  # @UnusedVariable
                csa_res_set = ResidueSet(name="csa", type="catalitic", residues=csas)
                strdoc.residue_sets.append(csa_res_set)


def update_quaternary(pdbUtils):
    '''
    Example – Author and computed assembly predictions agree
    REMARK 350 BIOMOLECULE: 1
    REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: DODECAMERIC
    REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE: DODECAMERIC
    
    Example – Author and computed assembly predictions differ
    REMARK 350 BIOMOLECULE: 1
REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: HEXAMERIC
REMARK 350 APPLY THE FOLLOWING TO CHAINS: A, B, C, D, E, F 
REMARK 350 BIOMOLECULE: 2
REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: HEXAMERIC
REMARK 350 APPLY THE FOLLOWING TO CHAINS: G, H, I, J, K, L
REMARK 350 BIOMOLECULE: 3
REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE: DODECAMERIC
REMARK 350 SOFTWARE USED: PISA
REMARK 350 TOTAL BURIED SURFACE AREA: 2990 ANGSTROM**2
REMARK 350 SURFACE AREA OF THE COMPLEX: 9330 ANGSTROM**2
REMARK 350 CHANGE IN SOLVENT FREE ENERGY: -40.0 KCAL/MOL
REMARK 350 APPLY THE FOLLOWING TO CHAINS: A, B, C, D, E, F, G, H, I,
REMARK 350 AND CHAINS: J, K, L
    '''
    total = ExperimentalStructure.objects().count()

    for strdoc in tqdm(ExperimentalStructure.objects().no_cache(), total=total):
        if not strdoc.quaternary:
            try:

                with open(pdbUtils.pdb_path(strdoc.name)) as h:
                    data = [l for l in h.readlines() if l.startswith("REMARK 350")]
                biomolecules_index = [i for i, l in enumerate(data) if "BIOMOLECULE:" in l] + [None]
                biomolecules = []
                for s, e in zip(biomolecules_index[0::2], biomolecules_index[1::2]):
                    biomolecule = data[s].split(":")[1].strip()
                    author = [l for l in data[s:e] if "AUTHOR DETERMINED BIOLOGICAL UNIT" in l]
                    if author:
                        author = author[0].split(":")[1].strip()
                    program = [l for l in data[s:e] if " SOFTWARE DETERMINED QUATERNARY STRUCTURE" in l]
                    if program:
                        program = program[0].split(":")[1].strip()
                    biomolecules.append((biomolecule, author, program))
                quaternaty = ""
                for bm in biomolecules:
                    quaternaty = "- Biomolecule " + str(bm[0]) + ": "
                    if (bm[1] or bm[2]) and (bm[1] == bm[2]):
                        quaternaty += ": " + bm[1]
                    elif bm[1]:
                        quaternaty += bm[1]
                    elif bm[2]:
                        quaternaty += bm[2]
                if len(biomolecules) == 1:
                    quaternaty = quaternaty.replace("- Biomolecule " + str(bm[0]) + ": ", "")
                strdoc.quaternary = quaternaty
                strdoc.save()
            except IndexError:
                _log.debug("no se puede parsear %s" % strdoc.name)
            except FileNotFoundError :
                _log.debug(f"{strdoc.name} could not be found")


def important_pfam(seqs_from_pdb_hmm):
    for query in tqdm(bpsio.parse(seqs_from_pdb_hmm, 'hmmer3-text')):
        try:
            pdb, chain, start, end = query.id.split("_")  # @UnusedVariable
            if ExperimentalStructure.objects(name=pdb,residue_sets__name="important_pfam").count():
                continue

            strdoc = ExperimentalStructure.objects(name=pdb).get()

            if not strdoc.residue_set("important_pfam"):
                important_rs = ResidueSet(name="important_pfam")
                domain_rs = None
                for hit in query:
                    if len(hit):
                        hsp = hit[0]
                        domain_rs = ResidueSet(name=hit.id)
                        i = 0
                        for x in str(hsp.aln[1].seq):
                            residue = chain + "_" + str(i + int(start))
                            if x == x.upper():
                                important_rs.residues.append(residue)
                            i = i + 1
                            domain_rs.residues.append(residue)
                        if domain_rs:
                            strdoc.residue_sets.append(domain_rs)
                strdoc.residue_sets.append(important_rs)
                strdoc.save()
        except DoesNotExist:
            pass


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

    parser.add_argument( "--csa", default='/data/databases/csa/csa.txt')
    parser.add_argument( "--hmm", default='/data/databases/pdb/pdb_seq_res.hmm')
    parser.add_argument( "--pdbs", default='/data/databases/pdb/')
    parser.add_argument( "--distances", default='/data/databases/pdb/processed/distances.tbl')


    args = parser.parse_args()


    #         pdbs = PDBs()
    #         pdbs.update('/data/pdb/divided/')

    BioMongoDB(args.database_name) #args.db_host

    # update_quaternary()
    #         # clusters cd hit
    #         update_clusters()
    #
    # residues near ligands --> metal drug/cofactor

    if not os.path.exists(args.csa):
        sys.stderr.write("%s not found. Download it from %s" % (
            args.csa,
            "http://www.ebi.ac.uk/thornton-srv/databases/CSA/downloads/CSA_2_0_121113.txt"
        ))
        sys.exit(1)

    if not os.path.exists(args.pdbs):
        sys.stderr.write("%s not found. Specify where is pdbs/divided directory" % (
            args.pdbs
        ))
        sys.exit(1)
    if not os.path.exists(args.distances):
        sys.stderr.write("%s not found. Run extended_domain.py script to create it." % (
            args.distances
        ))
        sys.exit(1)


    pdbUtils = PDBs(pdb_dir=args.pdbs)
    print("Update Quaternary")
    update_quaternary(pdbUtils)
    print("Update CSA")
    update_csa(args.csa)
    print("Update CYS/TYR")
    free_cys_tyr(pdbUtils)


    print("Update Importan Pfam")
    important_pfam(args.hmm)
    print("Update Binding residues")
    update_binding_residues(args.distances)
    _log.info("update pdb properties finished!!")




if __name__ == "__main__":
    sys.exit(main())
