"""
BEFORE installing  ete3
sudo apt-get install python-pyqt5 python-pyqt5.qtopengl python-pyqt5.qtsvggit s



"""


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
        newick = scipy2Newick(node.get_left(), newick, node.dist, leaf_names)
        newick = scipy2Newick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick

