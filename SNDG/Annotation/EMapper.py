from collections import defaultdict


class EMapper:

    def __init__(self):
        self.data = {}
        self.fields = "query_name	seed_eggNOG_ortholog	seed_ortholog_evalue	seed_ortholog_score	best_tax_level	Preferred_name	GOs	EC	KEGG_ko	KEGG_Pathway	KEGG_Module	KEGG_Reaction	KEGG_rclass	BRITE	KEGG_TC	CAZyBiGG_Reaction".split()

    def read_file(self, file_path):
        with open(file_path) as h:
            for line in h:
                vec = line.strip().split("\t")
                if not line.startswith("#"):
                    q = vec[0]
                    self.data[q] = {f: vec[i-1] for i, f in enumerate(self.fields, 1)}

    def group_prop(self, groups, prop_name):
        """
        :param groups: dict with group_name as key and list of proteins as value
        :param prop_name: field/property to be grouped, for example "EC"
        """
        prop_data = defaultdict(lambda: [])

        for g, prots in groups.items():
            for p in prots:
                prop_data[g] += self.data[p][prop_name].split(",")

        return {k: set(v) for k, v in prop_data.items()}
