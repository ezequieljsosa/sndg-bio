
def scipy2Newick(node, newick, parentdist, leaf_names):
    """
    https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format

    from scipy.cluster import hierarchy
    tree = hierarchy.to_tree(Z,False)
    newick = getNewick(tree, "", tree.dist, leaf_names)

    :param node:
    :param newick:
    :param parentdist:
    :param leaf_names:
    :return:
    """
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick

