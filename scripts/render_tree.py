#!/usr/bin/env python

'''
Created on Aug 24, 2015

@author: eze
'''
import os

os.environ['QT_QPA_PLATFORM'] = 'offscreen'
import argparse
import logging
from tqdm import tqdm
import pandas as pd
import warnings

from Bio import BiopythonWarning, BiopythonParserWarning

warnings.simplefilter('ignore', BiopythonWarning)
warnings.simplefilter('ignore', BiopythonParserWarning)

from SNDG.Comparative.Tree import TreeUtils

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", required=True, help="tree")
    parser.add_argument("--image", required=True, help="output image")
    parser.add_argument("--outgroup", default=None, help="coma separated outgroup nodes")
    parser.add_argument("--prune", default=None, help="coma separated nodes whose common ancestor will define the root")
    parser.add_argument("--metadata", action=None, help="metadata csv. first column must be the node name")
    parser.add_argument("--order", default=None, help="coma separated order of the metadata columns in the image")

    args = parser.parse_args()

    assert os.path.exists(args.tree), "%s does not exists" % args.tree

    if args.metadata:
        assert os.path.exists(args.metadata)

    tu = TreeUtils()
    tu.tree_path = args.tree
    tu.load_tree()
    metadata = None
    if args.metadata:
        data = {}
        metadata = pd.read_csv(args.metadata)
        node_id = metadata.columns[0]
        metadata = pd.read_csv(args.metadata, dtype={node_id: str})
        metadata = metadata.fillna("")

        for c in metadata.columns:
            metadata[c] = metadata[c].astype(str)
        for _, r in metadata.iterrows():
            data[r[node_id].split(".variant")[0]] = {c: r[c] for c in metadata.columns}

        metadata = data

    tu.render_tree(args.image, args.outgroup.split(",") if args.outgroup else None,
                   args.prune.split(",") if args.prune else None, metadata,
                   args.order.split(",") if args.order else None)
