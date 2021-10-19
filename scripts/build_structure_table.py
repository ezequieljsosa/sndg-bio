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
import ast
import logging
import os
import sys
from glob import glob
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import Bio.SeqIO as bpio
import Bio.SearchIO as bpsio
import pandas as pd
import subprocess as sp 
from SNDG import init_log

init_log()
_log = logging.getLogger(__name__)

__all__ = []
__version__ = 0.1
__date__ = '2020-09-18'
__updated__ = '2020-09-18'

from tqdm import tqdm


def main(argv=None):  # IGNORE:C0111
    from SNDG import mkdir
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

    parser.add_argument("-m", "--models_dir", required=True)
    # parser.add_argument("-a", "--alns_dir", required=True)
    parser.add_argument("-p", "--props_output", default="./structures.csv")
    parser.add_argument("-o", "--pdbs_output", default="./structures")

    parser.add_argument('-V', '--version', action='version', version=program_version_message)

    args = parser.parse_args()
    alns = list(glob(f"{args.models_dir}/*/*/*.ali"))
    data = []
    mkdir(args.pdbs_output)
    pbar = tqdm(alns)
    for aln in pbar:
        pbar.set_description(aln)
        q, h = bpio.parse(aln, "pir")
        # hsp = [hit for hit in bpsio.read(f"{args.alns_dir}/{q.id.replace('P1;', '')}.xml", "blast-xml") if
        #        hit.id.replace('P1;', '') == h.id.replace('P1;', '').replace("_", "")][0][0]
        pdbs = glob(aln.replace("seq.ali", "") + "*.pdb")
        if not pdbs:
            continue
        pdb = pdbs[0]
        model = pdb.split("/")[-2]
        with open("process.log","w") as h2:
            sp.call(f"ln {pdb} {args.pdbs_output}/{model}.pdb",shell =True,stderr=h2)
        with open(pdb + ".qmeans.json") as hqmean:
            xx = hqmean.read().replace("'", '"')
            try:
                qmean = ast.literal_eval(xx)
            except:
                sys.stderr.write(f"error reading qmeam from {pdb}\n")
                continue
        r = {
            "model":  q.id.replace("P1;", "") + "_" + h.id.replace("P1;", ""),
            "prot": q.id.replace("P1;", ""),
            "template": h.id.replace("P1;", ""),
            "qaln": str(q.seq),
            "haln": str(h.seq),
            "qstart": int(q.description.split(":")[2]),
            "qend": int(q.description.split(":")[2]) + int(q.description.split(":")[4][1:]),
            "hstartres": int(h.description.split(":")[2]),
            "hendres": int(h.description.split(":")[2]) + int(h.description.split(":")[4][1:]),
            # "hstart": hsp.hit_start,
            # "hend": hsp.hit_end,
            "qmean": qmean["qmean"],
            "zqmean": qmean["zqmean"],

        }
        data.append(r)

    pd.DataFrame(data).to_csv(args.props_output)


if __name__ == "__main__":
    sys.exit(main())
