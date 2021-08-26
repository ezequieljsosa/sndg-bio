from SNDG.BioMongo.Process.Importer import _common_annotations
import sys

import os
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

os.environ["COMPOUND_TYPES_PATH"] = os.getenv('COMPOUND_TYPES_PATH', "/target/data/compound_type.csv")
from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB

if __name__ == "__main__":
    argv = sys.argv

    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose", dest="verbose", action="count",
                        help="set verbosity level [default: %(default)s]")

    parser.add_argument("-host", "--db_host", default='127.0.0.1')
    parser.add_argument("--port", default=27017)
    parser.add_argument("-db", "--db_name", default='tdr')
    parser.add_argument("--pdbs_path", required=True)
    parser.add_argument("--organism_name", required=True)
    parser.add_argument("--remove_tmp", action='store_true')
    parser.add_argument("--cpu", default=4)
    parser.add_argument("--tmp_dir", default="./annotation/")

    args = parser.parse_args()

    mdb = BioMongoDB(args.db_name, port=args.port, host=args.db_host)
    _common_annotations(args.organism_name, args.tmp_dir, args.cpu, args.remove_tmp, True, False, None, args.pdbs_path)
