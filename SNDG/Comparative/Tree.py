"""
BEFORE installing  ete3
pip install pyqt5 pyqt5.qtopengl pyqt5.qtsvg


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


from ete3 import Tree, TreeStyle, TextFace
import vcf
from _collections import defaultdict
from tqdm import tqdm


def andDiff(d, conjuntos):
    x = set()
    y = d[conjuntos[0]]
    for c in conjuntos:
        x = x | d[c]
        y = y & d[c]
    return (len(y), len(x - y), [len(c - y) for c in conjuntos])


class TreeUtils(object):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.vcfs = []
        self.gvcf = "/tmp/samples.gvcf"
        self.tree_path = "/tmp/samples.newick"
        self.csv_shared_snps = "/tmp/shared_count.txt"
        self.sample_variants = defaultdict(lambda: [])
        self.shared_snps = defaultdict(lambda: 0)

    def load_tree(self):
        self.tree = Tree(self.tree_path)

    def tree_style_with_data(self,data={},order=None,force_topology = False):
        """
        newick: text or file
        render_in: default %%inline for notebooks
        """

        ts = TreeStyle()
        if data:
            if not order:

                first = data.keys()
                first = list(first)[0]
                order = data[first].keys()
            ts.show_leaf_name = True
            ts.draw_guiding_lines = True
            ts.force_topology = force_topology



        for i,x in enumerate(order):
            tf = TextFace(x)
            tf.margin_left = 5
            ts.aligned_header.add_face(tf, column=i)
        if data:
            for leaf in self.tree.get_leaves():
                for i,col in enumerate(order):
                    tf = TextFace(data[leaf.name][col])
                    tf.margin_left = 5
                    leaf.add_face(tf, column=i, position = "aligned")

        return ts

    def prune_from_common_ancestor(self,leafs_for_ca):
        ca = self.tree.get_common_ancestor(leafs_for_ca)
        self.tree.prune(ca)

    def render_tree(self,path, outgroup=None,prune=False,metadata=None,order=None):
        assert self.tree
        self.fix_gvcf_tree(outgroup)
        ts = self.tree_style_with_data(metadata,order) if metadata else None

        if prune:
            self.prune_from_common_ancestor(prune)
        self.tree.render(path,tree_style=ts)

    def fix_gvcf_tree(self, outgroup=None):

        for l in self.tree:
            l.name = l.name.split(".variant")[0]
        if outgroup:
            ancestor = self.tree.get_common_ancestor(*outgroup)
            if ancestor != self.tree:
                self.tree.set_outgroup(ancestor)

        for idx, node in enumerate(self.tree.traverse("levelorder")):

            if not node.name:
                node.name = "inner_" + str(idx)  # + ( ("_" + node.shared) if hasattr(node,"shared") else "")
        self.tree.name = "root"

    def load_gvcf_allele_dict(self, sample_name_fn=lambda x: x):
        with open(self.gvcf) as h:
            variantes = list(tqdm(vcf.VCFReader(h)))
        for v in tqdm(variantes):
            for sample in v.samples:
                if sample.data.GT != ".":
                    self.sample_variants[sample_name_fn(sample.sample)].append(
                        str(v.POS) + ":" + str(v.ALT[int(sample.data.GT) - 1]))
        self.sample_variants = dict(self.sample_variants)
        for x in self.sample_variants:
            self.sample_variants[x] = set(self.sample_variants[x])

    def export_exclusive_SNPs_count(self):
        with open(csv_shared_snps, "w") as h:
            for node in self.tree.traverse():
                count = self.shared_snps[node.name]
                h.write(node.name + " " + count + "\n")

    def init_shared_and_exclusive_SNPs(self, allow_missing=False):

        for node in self.tree.traverse():
            leafs = [x.name.replace(" ", "_") for x in (self.tree & node.name) if
                     x.name.replace(" ", "_") in self.sample_variants]
            if leafs:
                assert allow_missing or (len(self.tree) == len(leafs)), \
                    "gvcf samples and tree leafs  do not match: %i, %i" % (len(self.tree), len(leafs))
                intersection_set = set(self.sample_variants[leafs[0]])
                for descendant in leafs:
                    intersection_set = intersection_set & set(self.sample_variants[descendant])
                for descendant, alleles in self.sample_variants.items():
                    if descendant not in leafs:
                        intersection_set = intersection_set - alleles
                count = str(len(intersection_set))
                if node.name not in leafs:
                    self.shared_snps[node.name] = count
                    node.add_features(shared=count)


if __name__ == '__main__':
    treeUtils = TreeUtils()
    # treeUtils.gvcf = "/data/projects/23staphylo/processed/core-mapping-n315/strains.gvcf"
    # treeUtils.tree_path = "/data/projects/23staphylo/processed/core-tree-n315/vk-upgma.newick"
    # treeUtils.load_tree()
    #
    # treeUtils.load_gvcf_allele_dict()
    # treeUtils.init_shared_and_exclusive_SNPs(allow_missing=True)
    # treeUtils.fix_gvcf_tree(["1493", "1764", "0484"])
    #
    # treeUtils.tree.write(format=1,
    #                            outfile="/data/organismos/SaureusN315/projects/5ad8958af0c51a8f8d44cc95/vcfkit.newick")
    # print "OK"

    treeUtils = TreeUtils()
    treeUtils.tree_path = "./data/processed/diff_alelos/RAxML_bestTree.TEST"
    treeUtils.load_tree()
    treeUtils.fix_gvcf_tree()
    treeUtils.tree.set_outgroup("LP6005441-DNA_A08")
    treeUtils.fix_gvcf_tree()
    treeUtils.gvcf = "./data/processed/diff_alelos/var.vcf"
    treeUtils.load_gvcf_allele_dict()
    for k,v in treeUtils.sample_variants.items():
        treeUtils.sample_variants[k.split(".variant")[0]] = v
    treeUtils.init_shared_and_exclusive_SNPs(allow_missing=True)


    #     for x in treeUtils.tree:
    #         x.name = x.name.strip().split(".")[0]
    #     treeUtils.tree.write(format=1, outfile="/data/projects/Staphylococcus/annotation/tree/tree_vcf2.newick")
    #
    #     treeUtils.load_gvcf_allele_dict(sample_name_fn=lambda x:x.split(".")[0])
    # treeUtils.csv_shared_snps = "/media/eze/Data/data/projects/Staphylococcus/analysis/core_vc/high_qual/count_exclusive_SNPs.txt"
    # treeUtils.count_exclusive_SNPs(allow_missing=True,inner_prefix="")

    # treeUtils.tree.set_outgroup( self.tree&"NZ_LN854556" )
#     treeUtils.tree.show()
#     treeUtils.tree.write(format=1, outfile="/data/projects/Staphylococcus/annotation/tree/outgrouped2.newick")
