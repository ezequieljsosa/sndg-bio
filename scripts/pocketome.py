#!/usr/bin/env python

'''
Created on Aug 24, 2015

@author: eze
'''

import warnings

from Bio import BiopythonWarning, BiopythonParserWarning

warnings.simplefilter('ignore', BiopythonWarning)
warnings.simplefilter('ignore', BiopythonParserWarning)

import os
import argparse
import logging
from tqdm import tqdm
import glob

from SNDG import init_log
from SNDG.Structure.FPocket import FPocket

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--models_directory", required=True)
    parser.add_argument("-l", "--log_path", default=None)

    args = parser.parse_args()

    assert os.path.exists(args.models_directory), "%s does not exists" % args.models_directory

    if not args.log_path:
        args.log_path = args.models_directory + "/pocketome.log"
    init_log(args.log_path, logging.INFO)
    _log = logging.getLogger("pocketome")

    with tqdm(glob.glob(args.models_directory + "/*.pdb")) as pbar:
        for pdb in pbar:
            pbar.set_description(pdb)
            try:
                pocket_data = pdb.replace(".pdb","") + ".json"
                if not os.path.exists(pocket_data):
                    fpo = FPocket(pdb)
                    result = fpo.hunt_pockets()
                    result.save(pocket_data)
                    result.delete_dir()
            except Exception as e:
                _log.warn(e)
