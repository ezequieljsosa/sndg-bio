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
import sys
import logging
from tqdm import tqdm
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

import pymongo
from SNDG import init_log
from SNDG.BioMongo.Model.Protein import Protein
from SNDG.BioMongo.Process.StructureAnotator import StructureAnotator
from SNDG.BioMongo.Process.BioMongoDB import  BioMongoDB
from mongoengine.errors import DoesNotExist

from SNDG.BioMongo.Model.ResidueAln import ResidueAln

init_log()
_log = logging.getLogger(__name__)

__all__ = []
__version__ = 0.1
__date__ = '2015-09-18'
__updated__ = '2015-09-18'


def main(argv=None):  # IGNORE:C0111

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

    parser.add_argument("-n", "--name", required=True)
    parser.add_argument("-dir", "--structs_dir", required=True)
    parser.add_argument("-db_structure", "--db_structure",help="Mongo structure db", default='pdb')
    parser.add_argument("-db_genome", "--db_genome",help="Mongo proteins db", default='xomeq')
    parser.add_argument("-host", "--db_host", default='127.0.0.1')
    parser.add_argument('-V', '--version', action='version', version=program_version_message)

    args = parser.parse_args()

    BioMongoDB(args.db_genome)
    db = pymongo.MongoClient(args.db_host)[args.db_structure]

    sa = StructureAnotator(args.structs_dir + "/")
    total = sa.total(db, args.name, {})

    with tqdm(sa.iterator(db, args.name, {}), total=total) as pbar:
        for model in pbar:
            pbar.set_description(model.name)

            template = model.templates[0]

            try:
                protein = Protein.objects(organism=args.name, alias=template.aln_query.name).get()
            except DoesNotExist:
                _log.warn(template.aln_query.name + " does not exists")
            sa.annotate_model(model, protein.domains())
            model.save()


if __name__ == "__main__":
    sys.exit(main())
