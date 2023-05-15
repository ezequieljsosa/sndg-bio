#!/usr/bin/env python3

'''
Created on Aug 24, 2015

@author: eze
'''
import os

os.environ['QT_QPA_PLATFORM'] = 'offscreen'
import argparse
import pandas as pd

from ete3 import Tree, TreeStyle, TextFace, PhyloTree

if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--tree", required=True, help="tree in newick format")
    parser.add_argument("--image", required=True, help="output image")
    parser.add_argument("--metadata", action=None, help="metadata csv. first column must be the node name")
    parser.add_argument("--order", default=None, help="coma separated order of the metadata columns in the image")

    parser.add_argument("--outgroup", required=False, help="name of the node or list of nodes separated by coma "
                                                           "to select a common ancestor (internal node)")
    parser.add_argument("--remove", required=False,
                        help="nodes to remove separated by comma, this is done AFTER the outgroup is selected")
    parser.add_argument("--swap", required=False, nargs="+",
                        help="""Rotates the children of one or more internal nodes. 
                        Each internal node is described by its named children separated by comma,   
                        and each internal node is separated by space.
                         --swap  A,B  C,D   --> swaps the children of the common ancestor of A and B , 
                                                and the same for C and D
                          
                              /-A 
                           /-| 
                          |   \-B 
                        --| 
                          |   /-C 
                           \-| 
                             |   /-D 
                              \-| 
                                 \-H 
                        
                        After --swap  A,B  C,D               
                             /-B
                           /-|
                          |   \-A
                        --|
                          |      /-D
                          |   /-|
                           \-|   \-H
                             |
                              \-C
                               """)

    parser.add_argument("--force_topology", action='store_true')
    parser.add_argument("--show_length", action='store_true')
    parser.add_argument("--show_support", action='store_true')
    parser.add_argument("--show_scale", action='store_true')

    parser.add_argument("--circular", action='store_true')

    args = parser.parse_args()

    assert os.path.exists(args.tree), "%s does not exists" % args.tree

    if args.metadata:
        assert os.path.exists(args.metadata)

    tree = Tree(args.tree)

    data = {}
    order = []
    metadata = None
    if args.metadata:

        metadata = pd.read_csv(args.metadata)
        node_id = metadata.columns[0]
        metadata = pd.read_csv(args.metadata, dtype={node_id: str})
        metadata = metadata.fillna("")

        for c in metadata.columns:
            metadata[c] = metadata[c].astype(str)
        for _, r in metadata.iterrows():
            data[r[node_id].split(".variant")[0]] = {c: r[c] for c in metadata.columns}

        metadata = data

    # tu.render_tree(args.image, args.outgroup.split(",") if args.outgroup else None,
    #               args.prune.split(",") if args.prune else None, metadata,
    #               args.order.split(",") if args.order else None)

    ts = TreeStyle()

    if args.circular:
        ts.mode = "c"
    ts.force_topology = args.force_topology
    ts.show_branch_length = args.show_length
    ts.show_branch_support = args.show_support
    ts.show_scale = args.show_scale

    if args.outgroup:
        outgroup = tree.get_common_ancestor(*args.outgroup.split(",")) if "," in args.outgroup else args.outgroup
        tree.set_outgroup(outgroup)
    if args.swap:
        for s in args.swap:
            tree.get_common_ancestor(*s.split(",")).swap_children()

    if data:
        if not args.order:
            # order = list(data.keys())
            first = data.keys()
            first = list(first)[0]
            order = data[first].keys()
        else:
            order = args.order.split(",")
        ts.show_leaf_name = True
        ts.draw_guiding_lines = True


    for i, x in enumerate(order):
        tf = TextFace(x)
        tf.margin_left = 5
        ts.aligned_header.add_face(tf, column=i)
    if data:
        for leaf in tree.get_leaves():
            for i, col in enumerate(order):
                tf = TextFace(data[leaf.name][col])
                tf.margin_left = 5
                leaf.add_face(tf, column=i, position="aligned")

    tree.render(args.image, tree_style=ts)

    #  render_tree.py --tree pepe.newick --image pepe1.png
