"""
blastp -num_threads 10 -query proteinas.faa -db VFDB_setB_pro.fas.gz -outfmt 6 -qcov_hsp_perc 75 | sed 's|100.000|100|' | awk '$3>80'

"""
if __name__ == '__main__':
    import argparse
    import os
    import sys
    import pandas as pd
    import Bio.SeqIO as bpio
    from BCBio import GFF

    parser = argparse.ArgumentParser(description='Utils over genebank file')

    parser.add_argument('prokka_gff')
    parser.add_argument('blast_vf')
    parser.add_argument('rgi')

    # cmd.add_argument('output_gb', nargs='?', default=sys.stdout)

    args = parser.parse_args()

    lt_map = {}

    for contig in GFF.parse(args.prokka_gff):
        for f in contig.features:
            for sf in f.sub_features:
                if sf.type == "CDS" and "locus_tag" in sf.qualifiers:
                    lt_map[sf.qualifiers["locus_tag"][0]] = sf

    df_vf = pd.read_table(args.blast_vf, sep="\t",
                          names="qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore".split())
    for lt, df_vfs in df_vf[["qseqid", "sseqid","pident"]].groupby("qseqid"):


        r = lt_map[lt]
        prod = r.qualifiers.get("product", r.qualifiers["locus_tag"])[0]
        gene = prod.split("gene=")[1].split("]")[0] if "gene=" in prod else lt
        print(f'{lt} {gene} {" ".join([f"{x.sseqid}|{x.pident}" for x in df_vfs.itertuples()])}')
