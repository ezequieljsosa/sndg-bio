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
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import logging
import os, sys

from mongoengine.errors import DoesNotExist

from Bia import init_log
from Bia.Model import BioProperties
from Bia.Model.collections import Cluster

from Bia.Programs.Cluster.CDHit import CDHit

from BiaStructure.Model.structure import ExperimentalStructure, ResidueSet
from tqdm import tqdm

from Bio.PDB.PDBParser import PDBParser
import Bio.SearchIO as bpsio
import pandas as pd
from BiaStructure.IO.pdb import PDBs

init_log()
_log = logging.getLogger(__name__)

from SNDG.Structure.CompoundTypes import compound_type, get_compound_type

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
                cristal_doc = ExperimentalStructure.objects(name=pdb).get()
                cristal_doc.clusters = [x for x in cristal_doc.clusters if x.type != "PDB_Segments_95"]
                if not cristal_doc.cluster(cluster_name):
                    cristal_doc.clusters.append(cluster)
                    cristal_doc.save()
            except DoesNotExist as ex:
                print str(ex)
                pass


def update_drug_binding_residues():
    str_num = 0
    df_drug_dist = pd.read_table("/data/databases/pdb/processed/distances.tbl", sep="\t")
    for pdb, df_drug_dist_pdb in df_drug_dist.groupby("pdb"):
        str_num = str_num + 1
        _log.debug("procesing binding: " + pdb + " : " + str(str_num))
        # df_drug_dist_pdb_near = df_drug_dist[(df_drug_dist.pdb == strdoc.name) & (df_drug_dist.distance <= 4)  ]

        try:

            #             if not strdoc.residue_set("drug_binding"):
            df_drug_dist_pdb_near = df_drug_dist_pdb[df_drug_dist_pdb.distance <= 3]
            if len(df_drug_dist_pdb_near):
                strdoc = ExperimentalStructure.objects(name=pdb).get()
                residue_list = list(set([row.chain + "_" + str(row.prot_res) for i, row in
                                         df_drug_dist_pdb_near.iterrows()]))  # @UnusedVariable
                rs = ResidueSet(name="drug_binding", residues=residue_list)
                strdoc.residue_sets.append(rs)
                strdoc.save()
        #             df_metal_dist_pdb = df_metal_dist[(df_metal_dist.pdb == strdoc.name) & (df_metal_dist.distance <= 4)  ]
        #             if df_drug_dist_pdb :
        #                 rs = ResidueSet("metal_binding",residues=list(set([row.chain + "_" +  str(row.prot_res) for i , row in df_metal_dist_pdb.iterrows() ])))
        #                 strdoc.residue_sets.append(rs)
        except DoesNotExist:
            pass


def update_binding_residues():
    str_num = 0
    df_binding_dist = pd.read_table("/data/databases/pdb/processed/distances2.tbl", sep="\t")
    df_binding_dist["comptype"] = map(lambda x: compound_type[x], df_binding_dist.drugname)

    groups = df_binding_dist.groupby("pdb")
    total = len(groups)
    for pdb, df_binding_dist_pdb in groups:
        str_num = str_num + 1
        _log.debug("procesing binding: " + pdb + " : " + str(str_num) + " / " + str(total))
        try:
            for r_comp_type in ['LIPID', 'METAL', 'NUCLEOTIDE', 'SUGAR', "DRUG", "COFACTOR"]:
                if r_comp_type != "COFACTOR":
                    comp_type = r_comp_type
                else:
                    comp_type = "DRUG"

                df_comp_dist_pdb_near = df_binding_dist_pdb[
                    (df_binding_dist_pdb.distance <= 3) & (comp_type == df_binding_dist_pdb.comptype)]
                if len(df_comp_dist_pdb_near):
                    strdoc = ExperimentalStructure.objects(name=pdb).get()
                    residue_list = list(set([row.chain + "_" + str(row.prot_res) for i, row in
                                             df_comp_dist_pdb_near.iterrows()]))  # @UnusedVariable
                    rs = ResidueSet(name=comp_type.lower() + "_binding", residues=residue_list, type="binding")
                    strdoc.residue_sets.append(rs)
                    strdoc.save()

        except DoesNotExist:
            pass


def free_cys_tyr():
    parser = PDBParser(PERMISSIVE=1, QUIET=1)
    str_num = 0
    for strdoc in ExperimentalStructure.objects():
        str_num = str_num + 1
        _log.debug("procesing free cys/tyr: " + strdoc.name + " : " + str(str_num))
        if not (strdoc.residue_set("free_cys") or strdoc.residue_set("free_tyr")):
            bp_pdb = list(parser.get_structure(strdoc.name, strdoc.file_path()))[0]
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


def update_csa():
    df_csa = pd.read_csv("/data/databases/csa/csa.txt")

    str_num = 0
    for strdoc in ExperimentalStructure.objects():
        str_num = str_num + 1
        _log.debug("procesing  csa " + strdoc.name + " : " + str(str_num))
        pdb = strdoc.name
        pdb_csa = df_csa[df_csa["PDB ID"] == pdb]
        if len(pdb_csa) > 0:
            csas = [(row["CHAIN ID"]) + "_" + str(row["RESIDUE NUMBER"]) for i, row in
                    pdb_csa.iterrows()]  # @UnusedVariable
            csa_res_set = ResidueSet(name="csa", type="catalitic", residues=csas)
            strdoc.residue_sets.append(csa_res_set)


def update_quaternary():
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
    pdbUtils = PDBs()
    total = ExperimentalStructure.objects().count()

    for str_num, strdoc in tqdm(ExperimentalStructure.objects().no_cache(), total=total):

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


def important_pfam():
    for query in bpsio.parse('/data/databases/pdb/seqs_from_pdb.hmm', 'hmmer3-text'):
        try:
            pdb, chain, start, end = query.id.split("_")  # @UnusedVariable
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

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s

  Created by user_name on %s.
  Copyright 2015 BIA. All rights reserved.

  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
''' % (program_shortdesc, str(__date__))

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument("-v", "--verbose", dest="verbose", action="count",
                            help="set verbosity level [default: %(default)s]")

        # parser.add_argument("-dir", "--structs_dir", required = True )
        parser.add_argument("-db", "--database_name", default='xomeq')
        parser.add_argument("-host", "--db_host", default='127.0.0.1')
        parser.add_argument('-V', '--version', action='version', version=program_version_message)

        # Process arguments
        args = parser.parse_args()
        verbose = args.verbose

        if verbose > 0:
            print("Verbose mode on")

        #         pdbs = PDBs()
        #         pdbs.update('/data/pdb/divided/')

        from mongoengine import connect
        connect(args.database_name, host=args.db_host)
        update_quaternary()
        #         # clusters cd hit
        #         update_clusters()
        #
        # residues near ligands --> metal drug/cofactor
        update_binding_residues()

        update_csa()

        important_pfam()

        free_cys_tyr()

        _log.info("update pdb properties finished!!")

    except Exception:

        import traceback
        print(traceback.format_exc())

        return 2


if __name__ == "__main__":
    sys.exit(main())
