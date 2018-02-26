"""

"""
import os
import subprocess as sp

import pandas as pd
from tqdm import tqdm
from collections import defaultdict

from SNDG import execute
from SNDG.Comparative.VcfSnpeffIO import VcfSnpeffIO


class VariantSet():
    """
/usr/local/bin/gatk
#!/bin/bash
java -jar /opt/GATK/GenomeAnalysisTK.jar $@
    """

    @staticmethod
    def create_gvcf(vcfs_path_list, gvcf_path, ref_path):
        """

        :param vcfs_path_list: list of paths of the vcf files
        :param gvcf_path: gvcf to be created
        :param ref_path: fasta from the reference genome
        :return:
        """
        cmd_template = """
        gatk -T CombineVariants    -R {ref} {vcfs}   -o {out}    -genotypeMergeOptions UNIQUIFY
        """
        vcfs = " ".join(["--variant " + x for x in vcfs_path_list])
        cmd = cmd_template.format(vcfs=vcfs, out=gvcf_path, ref=ref_path)
        execute(cmd)

    def __init__(self, path_gvcf, reference=None):
        assert os.path.exists(path_gvcf), path_gvcf + " does not exists"
        self.total_variants = int(sp.check_output('grep -vc "^#" ' + path_gvcf, shell=True))
        self.reference = reference
        self.gvcf = VcfSnpeffIO.parse(path_gvcf)

    def build_table(self):

        rows = []
        for variant, effects in tqdm(self.gvcf, total=self.total_variants):
            effect = effects[0]
            vresult = {
                "chrom": variant.CHROM,
                "gene": effect.geneid,
                "ref": variant.REF,
                "pos": variant.POS - 1,
                "impact": effect.impact,
                "type": "&".join(effect.annotation)}
            alternatives = []
            for sample in variant.samples:
                sample_name = sample.sample.split(".variant")[0]
                if sample.called:
                    vresult[sample_name] = str(variant.ALT[int(sample.data.GT) - 1])
                    if effect.aa_pos:
                        vresult["aa_pos"] = effect.aa_pos
                        vresult["aa_ref"] = effect.aa_ref
                        vresult["aa_alt"] = effect.aa_alt
                else:
                    vresult[sample_name] = variant.REF
                alternatives.append(vresult[sample_name])
            if len(set(alternatives)) > 1:
                rows.append(vresult)
                # df_result = df_result.append(vresult, ignore_index=True)

        return pd.DataFrame(rows)

    def dist_variants(self, df_variants, samples):
        dist = defaultdict(lambda: defaultdict(lambda: 0))
        for _, row in df_variants.iterrows():
            for i, s1 in enumerate(samples):
                for j, s2 in enumerate(samples):
                    if i > j:
                        if row[s1] != row[s2]:
                            dist[s1][s2] += 1
        return dist

    def diff_variants(self, df_variants, samples):
        # diff = defaultdict(lambda: defaultdict(lambda: []))
        df = df_variants
        idx = df.pos != df.pos
        for i in range(len(samples)):
            s1 = df[samples[i]]
            s2 = df[samples[(i + 1) % len(samples)]]
            idx2 = s1 != s2
            idx = idx | idx2

        df = df[idx]
        """
        for _, row in tqdm(df.iterrows(), total=len(df)):
            for i, s1 in enumerate(samples):
                for j, s2 in enumerate(samples):
                    if i > j:
                        if row[s1] != row[s2]:
                            diff[s1][s2].append((row.chrom, row.pos, row.gene, row["type"], row[s1], row[s2]))
        """
        return df[["chrom", "pos", "gene", "type"] + samples]


if __name__ == '__main__':
    import glob
    from SNDG import init_log

    init_log()
    vcfs = glob.glob("/home/eze/workspace/git/msmegmatis_mut/data/processed/variant_call/**/*.vcf")
    vcfs = [x for x in vcfs if "ann" not in x]
    VariantSet.create_gvcf(vcfs, "/tmp/pepe.gvcf",
                                  "/home/eze/workspace/git/msmegmatis_mut/data/external/GCF_000283295.1_ASM28329v1_genomic.fna")
    pepe = VariantSet("/tmp/pepe.ann.gvcf",
                      "/home/eze/workspace/git/msmegmatis_mut/data/external/GCF_000283295.1_ASM28329v1_genomic.fna")
    df = pepe.build_table()
    df.to_csv("/tmp/pepe.csv", columns=["pos", "gene", "type", "ref"] +
                                       ["MUT-11","MUT-12","MUT-7","WT1"] + ["aa_pos", "aa_ref", "aa_alt"])
    print pepe
