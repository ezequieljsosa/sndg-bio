from collections import Counter


def process_line(line, samples, categories):
    _, pos, _, ref, alts = line.split()[:5]
    pos = int(pos)
    alts = [ref] + alts.split(",")
    gts = [x[0:3].replace("|", "/") for x in line.strip().split()[9:]]

    def gt2allele(gt):

        alleles = list(set(gt.split("/")))
        if alleles[0] == ".":
            return "N"
        if len(alleles) == 1:

            if alleles[0] == "0":
                return ""  # alts[int(alleles[0])]
            else:
                return "x"
        else:
            return "*"  # alts[int(alleles[0])] + "/" + alts[int(alleles[1])]

    ann = line.split()[7].split("ANN=")[1].split(",")[0].split("|")
    _, vtype, impact, gene, lt, _, _, _, _, genepos, protpos = ann[:11]

    new_line = f"{pos},{ref},{vtype},{impact},{gene},{lt},{genepos[2:]},{protpos[2:]},"

    tot_vars = 0

    for cat, value_map in categories.items():

        for value, csamples in sorted(value_map.items(), key=lambda x: x[0]):
            invariant = []
            for s in csamples:
                if s in samples:
                    idx = samples.index(s)
                    gt = gts[idx]
                    if gt not in ["0/0", "./."]:
                        invariant.append(s)
                        tot_vars += 1

            variant_val = round(len(invariant) / len(csamples), 2)
            new_line += str(variant_val) + ","
    for idx, s in enumerate(samples):
        new_line += gt2allele(gts[idx]) + ","

    return new_line


if __name__ == '__main__':
    import argparse
    from tqdm import tqdm
    import pandas as pd
    import sys
    from collections import defaultdict

    parser = argparse.ArgumentParser(description='creates a CSV with the variant conservation')
    parser.add_argument('vcf', help='annotated vcf')
    parser.add_argument('csv', help='csv with categories. first column must be named "sample"')

    args = parser.parse_args()

    df_features = pd.read_csv(args.csv)
    df_features["sample"] = df_features["sample"].astype(str)
    categories = defaultdict(lambda: defaultdict(list))
    cols = "pos,ref,vtype,impact,gene,lt,genepos,protpos".split(",")

    for col in df_features.columns:
        if col != "sample":
            for _, r in df_features.iterrows():
                categories[col][str(r[col])].append(r["sample"])
    categories = dict(categories)

    for cat, value_map in categories.items():
        for value, _ in sorted(value_map.items(), key=lambda x: x[0]):
            cols.append(value)

    with open(args.vcf) as h:
        for line in h:
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    samples = [x.strip() for x in line.split()[9:]]
                    cols += samples
                    if set(df_features["sample"]) - set(samples):
                        sys.stderr.write("some samples are in the csv properties file and not on the vcf file\n")
                        sys.stderr.write( ",".join(set(df_features["sample"]) - set(samples)) + "\n")
                    if set(samples) - set(df_features["sample"]) :
                        sys.stderr.write("some samples are in the vcf file and not in the csv properties file\n")
                        sys.stderr.write( ",".join( set(samples) - set(df_features["sample"])) + "\n")

                    sys.stdout.write(",".join(cols) + "\n")
                continue
            break

        for line in tqdm(h, file=sys.stderr):
            try:
                newline = process_line(line, samples, categories)
                sys.stdout.write(newline + "\n")


            except:
                sys.stderr.write(line)
                raise
