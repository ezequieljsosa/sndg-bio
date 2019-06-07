#!/usr/bin/env python

'''
Created on Aug 24, 2015

@author: eze
'''
import os
import argparse
import multiprocessing
import logging
from tqdm import tqdm

import warnings

from Bio import BiopythonWarning, BiopythonParserWarning

warnings.simplefilter('ignore', BiopythonWarning)
warnings.simplefilter('ignore', BiopythonParserWarning)

from SNDG.Comparative.VariantSet import VariantSet

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True, help="vcf or csv input")
    parser.add_argument("-o", required=True, help="aln fasta output")
    parser.add_argument("-r", default=None, help="fasta reference")
    parser.add_argument("-d","--depth", default=30, help="min_depth")
    parser.add_argument("-rebuild", action='store_true', help="rebuild csv")
    parser.add_argument("-csv", default="/tmp/variants.csv", help="filtered csv")

    args = parser.parse_args()

    assert os.path.exists(args.i), "%s does not exists" % args.i
    assert os.path.exists(args.r), "%s does not exists" % args.r

    if args.rebuild and os.path.exists(args.csv):
        os.remove(args.csv)


    vs = VariantSet(args.i)
    if args.i.endswith("csv"):
        df = pd.read_csv(args.i)
    else:
        df = vs.build_table(min_depth=int(args.depth))
        df.to_csv(args.csv)
    vs.aln(args.o, args.csv, args.r)
