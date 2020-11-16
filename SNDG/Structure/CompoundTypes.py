main_compound_types = ['DRUG','LIPID', 'METAL', 'NUCLEOTIDE', 'SUGAR']
COMPOUND_TYPE = [
                    "?",  "MODIFIED", "SOLVENT", "ORGANIC",  "COFACTOR", "RESIDUE"] + main_compound_types

compound_type = {}
with open('/data/databases/pdb/manually_curated/compound_type.csv') as handle:
    compound_type = {x.replace('"',"").strip().split(",")[0]:x.replace('"',"").strip().split(",")[1] if x.replace('"',"").strip().split(",")[1] else "?" for x in handle.readlines()}


def get_compound_type(residue):
    return compound_type[residue.get_resname().strip()] if residue.get_resname().strip() in compound_type else "?"


