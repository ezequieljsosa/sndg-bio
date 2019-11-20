from collections import defaultdict
from BCBio import GFF

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
                    self.data[q] = {f: vec[i - 1] for i, f in enumerate(self.fields, 1)}

                    self.data[q]["desc"] = " ".join(vec[len(self.fields):])

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

    def add2gff(self, gff_in, gff_out):
        rs = []
        with open(gff_in) as h:
            for seqrecord in GFF.parse(h):
                rs.append(seqrecord)
                for f in seqrecord.features:
                    for sf in f.sub_features:
                        if sf in self.data:
                            if self.data["GOs"]:
                                sf.qualifiers["GO"] = self.data["GOs"].split(",")
                            if self.data["EC"]:
                                sf.qualifiers["EC_number"] = self.data["EC"].split(",")
                            if self.data["KEGG_Pathway"]:
                                sf.qualifiers["KEGG_Pathway"] = self.data["KEGG_Pathway"].split(",")
                            if self.data["KEGG_ko"]:
                                sf.qualifiers["KEGG_ko"] = self.data["KEGG_ko"].split(",")
                            if self.data["KEGG_Reaction"]:
                                sf.qualifiers["KEGG_Reaction"] = self.data["KEGG_Reaction"].split(",")
                            if self.data["BRITE"]:
                                sf.qualifiers["BRITE"] = self.data["BRITE"]
                            if self.data["Preferred_name"]:
                                sf.qualifiers["gene"] = self.data["Preferred_name"]
                            sf.qualifiers["seed_eggNOG_ortholog"] = self.data["seed_eggNOG_ortholog"]
                            sf.qualifiers["note"] = self.data["desc"].replace(";",",")

        with open(gff_out, "w") as h:
            GFF.write(rs, h)
