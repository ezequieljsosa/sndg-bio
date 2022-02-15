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
import json
import logging
import os
import sys
import glob
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from pymongo import MongoClient
import pandas as pd
from Bio.PDB.PDBParser import PDBParser
from Bio.SeqUtils import seq1

from SNDG import init_log
from SNDG.BioMongo.Model.SeqCollection import SeqCollection
from SNDG.BioMongo.Model.Protein import Protein
from SNDG.BioMongo.Model.Alignment import AlnLine
from SNDG.BioMongo.Model.ResidueAln import ResidueAln
from SNDG.BioMongo.Model.Structure import ModeledStructure, Chain, Molecule, StructureQuality, \
    ResidueSet
from SNDG.BioMongo.Process.StructureAnotator import StructureAnotator

init_log()
_log = logging.getLogger(__name__)

__all__ = []
__version__ = 0.1
__date__ = '2015-09-18'
__updated__ = '2015-09-18'

from tqdm import tqdm


from SNDG.Structure.CompoundTypes import get_compound_type

def main(argv=None):  # IGNORE:C0111

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

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

    parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose", dest="verbose", action="count",
                        help="set verbosity level [default: %(default)s]")

    parser.add_argument("-g", "--genome", required=True)

    parser.add_argument("-dir", "--structs_dir", required=True)
    parser.add_argument("-prop", "--properties", required=True)

    parser.add_argument("--pocket_start_idx", default=0)

    parser.add_argument("-pname", "--pipeline", help="Pipeline name", default="SNDG")

    parser.add_argument("-db_prots", "--genome_database_name", default='test_database')
    parser.add_argument("-db", "--database_name", default='pdb')
    parser.add_argument("-host", "--db_host", default='127.0.0.1')
    parser.add_argument('-V', '--version', action='version', version=program_version_message)

    args = parser.parse_args()

    if not os.path.exists(args.properties):
        _log.error("%s does not exists" % args.properties)
        sys.exit(1)
    if not os.path.exists(args.structs_dir):
        _log.error("%s does not exists" % args.structs_dir)
        sys.exit(2)

    from mongoengine import connect, register_connection
    connect(args.genome_database_name, host=args.db_host)
    register_connection(args.database_name, "pdb")
    db = MongoClient()[args.genome_database_name]
    collection = db.sequence_collection.find_one({"name": args.genome}, {"name": 1})
    assert collection, "genome %s not found" % args.genome

    df_structs = pd.read_csv(args.properties)

    parser = PDBParser(PERMISSIVE=1, QUIET=1)
    model_files = glob.glob(args.structs_dir + "/*.pdb")
    models_count = len(model_files)
    _log.info("%i models will be processed for %s that were generated with pipeline %s" % (
        models_count, args.genome, args.pipeline))

    with tqdm(model_files) as pbar:
        for model_file in pbar:

            pbar.set_description("processing %s" % model_file)

            model_path = model_file
            model_name = os.path.basename(model_file).replace(".pdb", '')
            model_data = df_structs[df_structs.model == model_name]
            if len(model_data) == 0:
                _log.warn(model_name + " has no entry in the properties file, it will not be loaded ")
                continue
            model_data = model_data.iloc[0]

            prot = list(Protein.objects(organism=args.genome, alias__iexact=model_data["prot"]))
            if len(prot) == 0:
                _log.warn("Not found: " + str(model_data["prot"]))
            try:
                strdoc = process_model(args.structs_dir, args.pipeline, collection["name"], collection["_id"],
                                       model_data, model_name, model_path, parser)
                strdoc.save()
            except Exception as ex:
                _log.error(ex)

    _log.info("Finished")


def process_model(structs_dir, pipeline, seq_col_name, seq_col_id, model_data, model_name, model_path, parser):
    structure = parser.get_structure(model_name, model_path)
    strdoc = ModeledStructure(name=model_name, seq_collection_id=seq_col_id, pipeline=pipeline,
                              organism=seq_col_name)
    for x in ["qmean", "dope", "zdope", "zqmean", "ga341"]:
        if x in model_data:
            try:
                squality = StructureQuality(name=x, value=float(model_data[x]))
                strdoc.qualities.append(squality)
            except Exception as ex:
                _log.warn("error parseando %s en %s: " % (x, model_name) + str(ex))
    model = structure[0]
    for chain in model:
        chaindoc = Chain(name=chain.id)
        strdoc.chains.append(chaindoc)

        aln_h = model_data["haln"]
        aln_q = model_data["qaln"]
        prot_start = int(float(model_data["qstart"]))
        prot_end = int(float(model_data["qend"]))

        pdb_start = int(float(model_data["hstartres"]))
        pdb_end = int(float(model_data["hendres"]))

        hit_start =  0 # int(float(model_data["hstart"]))
        hit_end = - 1 # int(float(model_data["hend"]))

        template_name = model_data["template"]

        aln_query = AlnLine(name=model_data["prot"], seq=aln_q, start=prot_start - 1, end=prot_end - 1)
        aln_hit = AlnLine(name=template_name, seq=aln_h, start=hit_start, end=hit_end)
        template = ResidueAln(aln_query=aln_query, aln_hit=aln_hit,
                              query_res_start=1,
                              query_res_end=len(chain),
                              hit_res_start=pdb_start,
                              hit_res_end=pdb_end)
        strdoc.templates.append(template)

        for residue in chain:
            molecule = Molecule(resid=residue.id[1],
                                chain=chain.id,
                                compound=seq1(residue.get_resname()),
                                compound_type=get_compound_type(residue))
            if molecule.compound_type == 'RESIDUE':
                chaindoc.residues.append(molecule)
            else:
                if not [x for x in strdoc.ligands if
                        (
                                x.compound_type == molecule.compound_type) and molecule.compound_type == 'SOLVENT']:
                    strdoc.ligands.append(molecule)

    pockets_json = structs_dir + "/" + model_name + ".pdb.json"

    if os.path.exists(pockets_json):
        rss = StructureAnotator.pocket_residue_set(pockets_json,list(model.get_atoms()))
        strdoc.pockets = rss

    pockets_json = structs_dir + "/" + model_name + ".json"

    if os.path.exists(pockets_json):
        rss = StructureAnotator.pocket_residue_set(pockets_json,list(model.get_atoms()))
        strdoc.pockets = rss
    return strdoc


if __name__ == "__main__":
    sys.exit(main())
