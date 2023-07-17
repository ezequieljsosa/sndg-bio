from collections import Counter


def process_line(line,samples,categories):

    _, pos, _, ref, alts = line.split()[:5]
    pos = int(pos)
    alts = [ref] + alts.split(",")
    gts = [x[0] for x in line.split()[9:]]
    gts = ["N" if x[0] == "." else alts[int(x[0])] for x in gts]

    ann =   line.split()[7].split("ANN=")[1].split(",")[0].split("|")
    _,vtype,impact,gene,lt,_,_,_,_,genepos,protpos = ann[:11]

    new_line = f"{pos},{ref},{vtype},{impact},{gene},{lt},{genepos[2:]},{protpos[2:]},"

    tot_vars = 0

    for cat,value_map in categories.items():

        for value,csamples in sorted(value_map.items(),key=lambda x:x[0]):
            invariant = []
            for s in csamples:
                if s in samples:
                    idx = samples.index(s)
                    gt = gts[idx]
                    if gt != ref:
                        invariant.append(s)
                        tot_vars += 1

            variant_val = round(len(invariant) / len(csamples),2)
            new_line += str(variant_val) + ","


    return new_line



if __name__ == '__main__':
    import argparse
    from tqdm import tqdm
    import pandas as pd
    import sys
    from collections import defaultdict

    parser = argparse.ArgumentParser(description='creates a CSV with the variant conservation')
    parser.add_argument('vcf',  help='annotated vcf')
    parser.add_argument('csv',  help='csv with categories. first column must be named "sample"')

    args = parser.parse_args()

    df_features = pd.read_csv(args.csv)
    categories = defaultdict(lambda : defaultdict(list))
    cols = "pos,ref,vtype,impact,gene,lt,genepos,protpos".split(",")





    for col in df_features.columns:
        if col != "sample":
            for _,r in df_features.iterrows():
                categories[col][r[col]].append(r["sample"])
    categories = dict(categories)

    for cat,value_map in categories.items():
        for value,_ in sorted(value_map.items(),key=lambda x:x[0]):
            cols.append(value)

    sys.stdout.write(",".join(cols) + "\n")
    with open(args.vcf) as h:
        for line in h:
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    samples = [x.strip() for x in line.split()[9:]]

                continue
            break

        for line in tqdm(h, file=sys.stderr):
            try:
                newline = process_line(line,samples,categories)
                sys.stdout.write(newline + "\n")


            except:
                sys.stderr.write(line)
                raise
