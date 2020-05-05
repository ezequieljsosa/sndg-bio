"""

"""
import os
import subprocess as sp
import sys
import pandas as pd
from tqdm import tqdm
from collections import defaultdict
import glob
from SNDG import execute
from SNDG.Comparative.VcfSnpeffIO import VcfSnpeffIO
import pysam
import Bio.SeqIO as bpio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import Counter
import traceback
import numpy as np
import json


def validate_variant_fn(sample_data, min_depth=30):
    has_depth = hasattr(sample_data, "AD") and sample_data.DP
    depth = False
    if hasattr(sample_data, "AD") and isinstance(sample_data.AD, list):

        if sum(sample_data.AD) > 0:
            frequent = (sample_data.AD[1]) * 1.0 / sum(sample_data.AD) >= 0.75
        else:
            frequent = False
        depth = sum(sample_data.AD) >= min_depth

    elif hasattr(sample_data, "AD"):
        depth = sample_data.AD >= min_depth
        frequent = True
    else:
        frequent = False

    return has_depth and depth and frequent


class Empty:
    pass


class VariantSet:
    def __init__(self, df, samples, variant_phasing=None):
        """
        df : pandas dataframe build by VariantSetUtils
        """
        self.df = df
        self.samples = samples
        self.variant_phasing = variant_phasing

    def sets(self,tqdm_fn=tqdm):
        sets = defaultdict(list)
        for _, r in tqdm_fn(self.df.iterrows(),total=len(self.df)):
            # alt = list(set([r[sample] for sample in self.samples]) - set(r.ref))[0]
            variant_id = "_".join([r.ref, str(r.pos)]) + "_"
            for sample in self.samples:
                if r[sample] != r.ref:
                    gta = VCFGenotypeAlleles.parse(r[sample] ,r.ref)
                    for allele in gta.alleles():
                        sets[sample].append(variant_id + allele)

        return {k: set(v) for k, v in sets.items()}

    def pos_analysis(self, low_cov=30, low_maf=0.75, include_col=None):
        report = []
        for _, row in self.df.iterrows():
            # pos = "_".join([row["chrom"], str(row["pos"]), row["ref"]])
            alts = [row[sample] for sample in self.samples
                    if row[sample] not in [row["ref"], "/".join([row["ref"], row["ref"]])]]
            complete = len(alts) == len(self.samples)

            with_low_depth = [(row[sample + "_AD_r"] + row[sample + "_AD_a"])
                              < low_cov for sample in self.samples]

            with_low_cov_alt = [(row[sample + "_AD_a"])
                                < low_cov for sample in self.samples
                                if row[sample] not in [row["ref"], "/".join([row["ref"], row["ref"]])]]

            with_low_cov_ref = [(row[sample + "_AD_r"] + row[sample + "_AD_a"])
                                < low_cov for sample in self.samples
                                if row[sample] not in [row["ref"], "/".join([row["ref"], row["ref"]])]]

            low_maf = any([(row[sample + "_AD_a"] * 1.0 /
                            (row[sample + "_AD_a"] + row[sample + "_AD_r"])
                            < low_maf) if (row[sample + "_AD_a"] + row[sample + "_AD_r"]) else True
                           for sample in self.samples
                           ])

            line = {
                "chrom": row["chrom"],
                "pos": row["pos"],
                "ref": row["ref"],
                "het": any([len(set(row[sample].split("/"))) == 2
                            for sample in self.samples if row["ref"] != row[sample]]),
                "low_maf": low_maf,
                "sample_with_alt": len(alts),
                "allele_warn": len(set(alts)) > 2,
                "complete": complete,
                "in_all": complete & (len(set(alts)) == 1),

                "low_depth_count": len([x for x in with_low_depth if x]),
                "low_cov_alt_count": len([x for x in with_low_cov_alt if x]),
                "low_cov_ref_count": len([x for x in with_low_cov_ref if x]),

                "low_depth": any(with_low_depth),
                "low_cov_alt": all(with_low_cov_alt),
                "low_cov_ref": all(with_low_cov_ref),

                "in_phase": (row["pos"] in self.variant_phasing) and (self.variant_phasing[row["pos"]] != row["pos"]),

            }

            report.append(line)
        report = pd.DataFrame(report)
        cols = self.df.columns

        if not include_col:
            cols = ["chrom", "pos", "ref", "indel"] + include_col + self.samples
        return pd.merge(self.df[cols], report, how='left', on=["chrom", "pos", "ref"])


class VCFGenotype():
    def __init__(self, gt):
        self.gt_str = gt
        self.gt = [x for x in self.gt_str.replace("|", "/").split("/")]
        self.ploidy = len(self.gt)

    def N(self):
        return "/".join("N" * self.ploidy)

    def isN(self):
        return self.N() == self.gt_str.replace(".","N").replace("|","/")

    def withAllele(self, ref, alts):
        return VCFGenotypeAlleles(self.gt_str, ref, alts)





class VCFGenotypeAlleles(VCFGenotype):
    def __init__(self, gt, ref, alts):
        super(VCFGenotypeAlleles, self).__init__(gt)
        self.ref = ref
        self.alts = [x for x in alts]
        assert self.ref not in alts
        all_alleles = self.alleles()
        self.allele_gt = [all_alleles[int(haplo)] if haplo != "." else "N" for haplo in self.gt]

    def __str__(self):
        return "/".join(self.allele_gt)

    def alleles(self):
        return [self.ref] + self.alts

    def alt_ploidy(self,allele=False):
        if allele:
            return len([x for x in self.allele_gt.split("/") if x == allele])
        else:
            return len([x for x in self.allele_gt.split("/") if x != self.ref])

    def is_ref(self):
        return len(set(self.alts) - set([self.ref])) == 0



    @staticmethod
    def parse(gt_with_allele_str,ref):
        alts = [x for x in gt_with_allele_str.split("/") if x != ref]
        alleles = [ref] + alts
        gt_str = "/".join([ str(alleles.index(x)) for x in gt_with_allele_str.split("/")])
        gt = VCFGenotypeAlleles(gt_str,ref,alts)
        return gt



class VariantSetUtils():
    """
bams_dict = {'1_20079_S5': '/home/andres/cepas/20079_S5_L001/dedup_reads.bam',
             '2_20651_S6': '/home/andres/cepas/20651_S6_L001/dedup_reads.bam',
            ...}

gvcf = VariantSet(path_vcf + 'combined.ann.gvcf',
                  path_vcf + 'h37.fna',bams_dict=bams_dict)

def select_default(variant, sample, alt, assigned_values):
    samples = [s.sample.split(".variant")[0] for s in variant.samples]
    sample_name = sample.sample.split(".variant")[0]
    idx = samples.index(sample_name)
    if idx - 1 == -1:
        return variant.REF
    else:
        return assigned_values[samples[idx - 1]]

gvcf.default_value_fn = select_default
df = gvcf.build_table()
    """

    @staticmethod
    def mark_with_vcf(df, bedfile, tag):
        df[tag] = [False for _ in range(len(df))]
        for line in open(bedfile):
            if line.startswith("#"): continue
            chrom, pos = line.split("\t")[:2]
            pos = int(pos)
            df.at[(df.chrom == chrom) & (df.pos == pos), tag] = True

    @staticmethod
    def correct_by_ad(vs, min_cov=20):
        df = vs.df.copy()
        for i, row in tqdm(df.iterrows(), total=len(df)):
            for sample in vs.samples:
                if ((row[sample] != row["ref"]) and
                        (row[sample + "_AD_r"] > row[sample + "_AD_a"]) and
                        (min_cov < row[sample + "_AD_a"])
                ):
                    df.loc[i, sample] = row["ref"]

        return VariantSet(df, vs.samples,
                          variant_phasing=vs.variant_phasing)

    @staticmethod
    def complete_vs_with_bamreadcount(df, sample, bamreadcountfile):
        """
        https://github.com/genome/bam-readcount
        NC_000962.3	1977	A	57	=:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	A:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	C:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	G:57:60.00:34.93:0.00:22:35:0.49:0.01:47.32:22:0.56:204.33:0.53T:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00	N:0:0.00:0.00:0.00:0:0:0.00:0.00:0.00:0:0.00:0.00:0.00
        """
        for line in tqdm(open(bamreadcountfile)):
            try:
                chrom, pos, base, deph, _, A, C, G, T = line.strip().split("\t")[:9]
                pos = int(pos)
            except:
                print(line)
                raise
            counts = {b: int(c.split(":")[1]) for b, c in zip(["A", "C", "G", "T"], [A, C, G, T])}

            qs = df[sample + "_AD_r"].isna() & (df.chrom == chrom) & (df.pos == pos)
            row = df[qs]

            if len(row):
                row = row.iloc[0]
                df.at[qs, sample + "_AD_r"] = counts[row.ref[0]]
                df.at[qs, sample + "_AD_a"] = counts[row[sample].split("/")[1][0]] if len(
                    row[sample].split("/")) > 1 else 0
        df[sample + "_AD_r"] = df[sample + "_AD_r"].fillna(0)
        df[sample + "_AD_a"] = df[sample + "_AD_a"].fillna(0)

    @staticmethod
    def create_gvcf(vcfs_folder, output_gvcf, ref_path):
        """

        :param vcfs_path_list: list of paths of the vcf files
        :param gvcf_path: gvcf to be created
        :param ref_path: fasta from the reference genome
        :return:
        """
        cmd_template = """
        docker run --rm -w /out {mount2} -v {vcfs_folder}:/out/vcfs/ -v {ref_folder}:/out/ref/ broadinstitute/gatk3:3.8-1 \
            java -jar /usr/GenomeAnalysisTK.jar -T CombineVariants    -R /out/ref/{ref_file} {vcfs} \
            -o {out_path}/{out_file} -genotypeMergeOptions UNIQUIFY
        """

        ref_folder = os.path.dirname(ref_path)
        ref_file = os.path.basename(ref_path)
        out_folder = os.path.dirname(output_gvcf)
        out_file = os.path.basename(output_gvcf) + ".bk"

        vcfs_path = "/out/vcfs/"
        if vcfs_folder == out_folder:
            mount2 = ""
            out_path = "/out/vcfs/"
        else:
            out_path = "/out/out/"
            mount2 = " -v {out_folder}:/out/out/ ".format(out_folder=out_folder)

        vcfs = " ".join(["--variant {vcfs_path}".format(vcfs_path=vcfs_path) + x for x in os.listdir(vcfs_folder) if
                         x.endswith(".vcf") or x.endswith(".vcf.gz")])
        cmd = cmd_template.format(vcfs=vcfs, out_folder=out_folder, out_file=out_file, mount2=mount2, out_path=out_path,
                                  ref_folder=ref_folder, ref_file=ref_file, vcfs_folder=vcfs_folder)
        print(cmd)
        execute(cmd)
        with open(out_folder + "/" + out_file) as h, open(output_gvcf, "w") as hw:
            for l in h:
                if l.startswith("#CHROM"):
                    vec = l.split("\t")
                    l = "\t".join(vec[:9] + [x.split(".variant")[0] for x in vec[9:]]) + "\n"

                hw.write(l)

    @staticmethod
    def haplotype_call(bam_path, output_gvcf, ref_path, ploidy=2, only_cmd=False):
        bam_folder = os.path.dirname(bam_path)
        bam_file = os.path.basename(bam_path)

        ref_folder = os.path.dirname(ref_path)
        ref_file = os.path.basename(ref_path)
        out_folder = os.path.dirname(output_gvcf)
        out_file = os.path.basename(output_gvcf)

        docker_bam_folder = "/out/bam/"
        if bam_folder == out_folder:
            mount2 = ""
            docker_out_path = docker_bam_folder
        else:
            docker_out_path = "/out/out/"
            mount2 = f" -v {out_folder}:{docker_out_path} "

        cmd = f"""docker run --rm -w /out {mount2} -v {bam_folder}:/out/bam/ -v {ref_folder}:/out/ref/ broadinstitute/gatk:4.1.0.0 \
        java -jar /gatk/gatk-package-4.1.0.0-local.jar  HaplotypeCaller -ERC GVCF \
        -R /out/ref/{ref_file} -ploidy {ploidy} \
        -I /out/bam/{bam_file} --output-mode EMIT_ALL_SITES \
        -O {docker_out_path}/{out_file}"""
        if only_cmd:
            return cmd
        else:
            execute(cmd)

    @staticmethod
    def combineGVCFs(vcfs_folder, output_gvcf, ref_path):
        """

        :param vcfs_path_list: list of paths of the vcf files
        :param gvcf_path: gvcf to be created
        :param ref_path: fasta from the reference genome
        :return:
        """
        cmd_template = """
        docker run --rm -w /out {mount2} -v {vcfs_folder}:/out/vcfs/ -v {ref_folder}:/out/ref/ broadinstitute/gatk:4.1.0.0 \
            java -jar /gatk/gatk-package-4.1.0.0-local.jar  CombineGVCFs    -R /out/ref/{ref_file} {vcfs} \
            -O {out_path}/{out_file} 
        """

        ref_folder = os.path.dirname(ref_path)
        ref_file = os.path.basename(ref_path)
        out_folder = os.path.dirname(output_gvcf)
        out_file = os.path.basename(output_gvcf) + ".bk"

        vcfs_path = "/out/vcfs/"
        if vcfs_folder == out_folder:
            mount2 = ""
            out_path = "/out/vcfs/"
        else:
            out_path = "/out/out/"
            mount2 = " -v {out_folder}:/out/out/ ".format(out_folder=out_folder)

        vcfs = " ".join(["--variant {vcfs_path}".format(vcfs_path=vcfs_path) + x for x in os.listdir(vcfs_folder) if
                         x.endswith(".vcf") or x.endswith(".vcf.gz")])
        cmd = cmd_template.format(vcfs=vcfs, out_folder=out_folder, out_file=out_file, mount2=mount2, out_path=out_path,
                                  ref_folder=ref_folder, ref_file=ref_file, vcfs_folder=vcfs_folder)
        print(cmd)
        execute(cmd)
        with open(out_folder + "/" + out_file) as h, open(output_gvcf, "w") as hw:
            for l in h:
                if l.startswith("#CHROM"):
                    vec = l.split("\t")
                    l = "\t".join(vec[:9] + [x.split(".variant")[0] for x in vec[9:]]) + "\n"

                hw.write(l)

    def __init__(self, path_gvcf, reference=None, bams_dict=None):
        assert os.path.exists(path_gvcf), path_gvcf + " does not exists"
        self.total_variants = int(sp.check_output('grep -vc "^#" ' + path_gvcf, shell=True))
        self.reference = reference
        self.bam = bams_dict
        if bams_dict:
            bams_dict2 = {}
            for sample, path in bams_dict.items():
                bams_dict2[sample] = pysam.AlignmentFile(path, "rb")
            self.bam = bams_dict2
        self.gvcf = VcfSnpeffIO.parse(path_gvcf)
        self.validate_variant_fn = validate_variant_fn
        self.default_value_fn = lambda variant, sample, alt, assigned_values: alt
        self.default_not_called_fn = lambda cov_ref, cov_alt: cov_ref < 30 or (
                (cov_ref * 1.0 / (cov_alt + cov_ref)) < 0.75)

    def build(self, variant_phasing=None, tqdm_fn=tqdm):

        variant_phasing_data = defaultdict(lambda: [])

        rows = []
        for variant, effects in tqdm_fn(self.gvcf, total=self.total_variants):
            if effects:
                effect = effects[0]
            else:
                effect = Empty()
                effect.geneid = ""
                effect.impact = ""
                effect.aa_pos = ""
                effect.aa_ref = ""
                effect.aa_alt = ""
                effect.hgvs_c = ""
                effect.hgvs_p = ""
                effect.annotation = []

            vresult = {
                "chrom": variant.CHROM,
                "gene": effect.geneid,
                "ref": variant.REF,
                "pos": variant.POS,
                "impact": effect.impact,
                "type": "&".join(effect.annotation),
                "qc_qual": variant.QUAL,
                "indel": (len(variant.REF) > 1 or any([len(str(x)) > 1 for x in variant.ALT]))
            }
            for qc_field in ["AC", "AF", "AN", "BaseQRankSum", "DP", "ExcessHet", "FS", "MLEAC", "MLEAF", "MQ",
                             "MQRankSum", "QD",
                             "MReadPosRankSum", "SOR"]:
                if qc_field in variant.INFO:
                    if isinstance(variant.INFO[qc_field], list):
                        vresult["qc_" + qc_field] = float(variant.INFO[qc_field][0])
                    else:
                        vresult["qc_" + qc_field] = float(variant.INFO[qc_field])

            alternatives = []
            for sample in variant.samples:
                sample_name = sample.sample.split(".variant")[0]
                gt = VCFGenotype(str(sample.data.GT))

                if sample.called:

                    if hasattr(sample.data, "PS") and sample.data.PS:  # should be PID
                        variant_phasing_data[int(sample.data.PS)].append(sample_name)

                    assert not gt.isN()

                    gt = gt.withAllele(variant.REF, [str(x) for x in variant.ALT])

                    alt = str(gt)

                    for qc_gt in ["AD", "GQ", "PL"]:
                        if hasattr(sample.data, qc_gt) and isinstance(getattr(sample.data, qc_gt), list):
                            vresult[sample_name + "_" + qc_gt + "_a"] = float(getattr(sample.data, qc_gt)[1])
                            vresult[sample_name + "_" + qc_gt + "_r"] = float(getattr(sample.data, qc_gt)[0])
                        elif hasattr(sample.data, qc_gt):
                            vresult[sample_name + "_" + qc_gt] = float(getattr(sample.data, qc_gt))
                        else:
                            vresult[sample_name + "_" + qc_gt] = np.NaN

                    vresult[sample_name] = alt

                    if (alt == vresult[sample_name]) and effect.aa_pos:
                        vresult["aa_pos"] = effect.aa_pos
                        vresult["aa_ref"] = effect.aa_ref
                        vresult["aa_alt"] = effect.aa_alt
                        vresult["gene_pos"] = effect.c_dna_pos
                elif gt.isN():
                    vresult[sample_name] = "N"
                else:
                    traceback.print_stack()
                    print(gt.N(),gt.isN(),gt)
                    raise Exception("Not called and invalid GT: " + str(sample.data.GT))

                alternatives.append(vresult[sample_name])

            rows.append(vresult)

        variant_phasing_data = dict(variant_phasing_data)

        if variant_phasing:
            with open(variant_phasing, "w") as h:
                json.dump(variant_phasing_data, h)
        return VariantSet(pd.DataFrame(rows),
                          [sample.sample.split(".variant")[0] for sample in variant.samples],
                          variant_phasing=variant_phasing_data)

    def dist_variants(self, df_variants, samples):
        dist = defaultdict(lambda: defaultdict(lambda: 0))
        for _, row in df_variants.iterrows():
            for i, s1 in enumerate(samples):
                for j, s2 in enumerate(samples):
                    if i > j:
                        if row[s1] != row[s2]:
                            dist[s1][s2] += 1
        return dist

    def diff_variants(self, df_variants, samples, columns=["chrom", "pos", "gene", "type"]):
        # diff = defaultdict(lambda: defaultdict(lambda: []))
        df = df_variants
        idx = df.pos != df.pos
        for i in range(len(samples)):
            s1 = df[samples[i]]
            s2 = df[samples[(i + 1) % len(samples)]]
            idx2 = s1 != s2
            idx = idx | idx2

        df = df[idx]

        return df[columns + samples]

    def aln(self, aln, txt="/tmp/samples.txt", ref=None):

        ref = str(bpio.read(ref, "fasta").seq)
        if not os.path.exists(txt):
            df = self.build_table()
            df.to_csv(txt)
        else:
            df = pd.read_csv(txt)
            df = df.fillna("")
        samples = [c for c in df.columns if c not in ["chrom", "pos", "ref", 'Unnamed: 0', "gene", "impact", "type", ]
                   and ("_ada" not in c) and ("_adr" not in c)]

        seqmap = {s: "" for s in samples}
        total = len(df)
        base_idx = 0
        for _, r in tqdm(df.sort_values("pos").iterrows(), total=total):
            pos_size = max([len(r[x]) for x in samples])
            for s in samples:
                subseq = ref[base_idx:r.pos] + r[s].ljust(pos_size, "-")
                seqmap[s] += subseq
            base_idx = r.pos + len(r.ref)

        for s in samples:
            seqmap[s] += ref[base_idx:]

        with open(aln, "w") as h:
            for k, v in seqmap.items():
                bpio.write(SeqRecord(id=k, name="", description="", seq=Seq(v)), h, "fasta")


if __name__ == '__main__':
    import glob
    from SNDG import init_log

    """
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
    
    
     import vcf
     ...: seqs = defaultdict(lambda:"")
     ...: alelos = []
     ...: #with open("./Full2mergeall_onlyvariants_maf.recode.vcf") as h:
     ...: #    variantes = list(vcf.VCFReader(h))
     ...: terminar = False
     ...: for v in tqdm(variantes):
     ...:         pos_size = [v.REF]
     ...:         alts = [v.REF] + [str(x) for x in v.ALT]
     ...:         alelos = []
     ...:         for sample in  v.samples:
     ...:             if (sample.data.GT != ".") and  (str(alts[int(sample.data.GT) ]) != "<CN0>"):
     ...:                 pos_size.append(  str(alts[int(sample.data.GT) ])  )
     ...:             if (sample.data.GT != "."):
     ...:                 alelos.append(str(alts[int(sample.data.GT) ]) )
     ...:         
     ...:                           
     ...:         if "<CN0>" in alelos:
     ...:             continue
     ...:         pos_size = max([len(x) for x in pos_size])
     ...:         for sample in v.samples:
     ...:             if sample.called:
     ...:                 alelo = str(alts[int(sample.data.GT) ])                
     ...:             elif sample.data.GT == ".":
     ...:                 alelo = "N"
     ...:             else:
     ...:                 alelo = v.REF            
     ...:             if alelo == "<CN0>":
     ...:                 alelo = ""
     ...:             alelo = alelo.ljust(pos_size, "-")
     ...:             seqs[sample.sample] += alelo
     ...:             if sample.sample in fasta:
     ...:                 if seqs[sample.sample] != str(fasta[sample.sample].seq)[: len(seqs[sample.sample])]:                
     ...:                     terminar = True
     ...:         if terminar:
     ...:             break
     ...: print seqs[sample.sample][-10:] 
     ...: print str(fasta[sample.sample].seq)[ len(seqs[sample.sample]) - 10: len(seqs[sample.sample])]

    
    """


    # strains = set(['0058', '1300', '0271', '0037', '1875', '3296', '1710', '1584',
    #                '1445', '0142', '1527', '0564', '3867NE', '1096', '0484', '1796',
    #                '3867NI', '1493', '0298', '1707', '3867INF', '1764', '0450'])
    #
    #
    # bams = {
    #     bam.split("/")[-2]:pysam.AlignmentFile(bam, "rb")
    #     for bam in glob("/mnt/data2/data/projects/23staphylo/data/steps/03-mappingn315/*/aln.bam")
    # }
    #
    #
    #
    # pepe = VariantSet("/mnt/data2/data/projects/23staphylo/data/steps/s07_joined_variant_call/output.ann.vcf",
    #                   "/mnt/data2/data/projects/23staphylo/external/n315/refN315.fasta",
    #                   bams_dict=bams)
    # df = pepe.build_table()

    def default_value_fn(variant, sample, alt, assigned_values):
        if sample.sample == "3867INF":
            s = [x for x in variant.samples if x.sample == "3867NI"][0]
            if s.called and validate_variant_fn(s.data):
                alt = str(variant.ALT[int(s.data.GT) - 1])
                return alt
            else:
                s = [x for x in variant.samples if x.sample == "3867NE"][0]
                if s.called and validate_variant_fn(s.data):
                    alt = str(variant.ALT[int(s.data.GT) - 1])
                    return alt
                return variant.REF
        elif sample.sample == "3867NE":
            s = [x for x in variant.samples if x.sample == "3867INF"][0]
            if s.called and validate_variant_fn(s.data):
                alt = str(variant.ALT[int(s.data.GT) - 1])
                return alt
            else:
                s = [x for x in variant.samples if x.sample == "3867NI"][0]
                if s.called and validate_variant_fn(s.data):
                    alt = str(variant.ALT[int(s.data.GT) - 1])
                    return alt
                return variant.REF
        elif sample.sample == "3867NI":
            s = [x for x in variant.samples if x.sample == "3867INF"][0]
            if s.called and validate_variant_fn(s.data):
                return str(variant.ALT[int(s.data.GT) - 1])
            else:
                s = [x for x in variant.samples if x.sample == "3867INF"][0]
                if s.called and validate_variant_fn(s.data):
                    return str(variant.ALT[int(s.data.GT) - 1])
                return variant.REF
        return alt


    # vs = VariantSet("/mnt/data2/data/projects/23staphylo_old/external/sordelli/variant_call/cohort.g.vcf")
    # # vs.default_not_called_fn = default_not_called_fn
    # vs.default_value_fn = default_value_fn
    #
    # df = vs.build_table()
    with open("/tmp/imputed.vcf", "w") as h:
        VariantSet.complete_vcf_with_bams("/home/eze/projects/ST5/processed/groups/arg_chi/all.vcf",
                                          "/home/eze/projects/ST5/processed/strains/*/alignment/mapped_reads.bam",
                                          stdout=h)
