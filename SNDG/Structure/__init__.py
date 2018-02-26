"""

depends

"""


def read_pdb_entries(pdb_entries_path):
    # entries_columns =  ["IDCODE", "HEADER", "ACCESSIONDATE", "COMPOUND", "SOURCE", "AUTHORS", "RESOLUTION", "EXPERIMENT"]
    # pdb_entries_df = pd.read_table(args.entries,skiprows=[0,1,2],sep='\t',names=entries_columns)
    def res_fn(y):
        try:
            res = float(y)
        except:
            res = 30
        return res


    entries = {x.split("\t")[0].lower(): res_fn(x.split("\t")[6]) for x in list(open(pdb_entries_path))[3:]}
    return entries




