"""

"""
import json
import os
import subprocess as sp
import sys
import traceback
from collections import defaultdict
from glob import glob
from itertools import groupby
import Bio.SeqIO as bpio
import numpy as np
import pandas as pd
import vcf
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from SNDG import execute
from SNDG.Comparative.VcfSnpeffIO import VcfSnpeffIO


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

    def sets(self, tqdm_fn=tqdm):
        sets = defaultdict(list)
        for _, r in tqdm_fn(self.df.iterrows(), total=len(self.df)):
            # alt = list(set([r[sample] for sample in self.samples]) - set(r.ref))[0]
            variant_id = "_".join([r.ref, str(r.pos)]) + "_"
            for sample in self.samples:
                if r[sample] != r.ref:
                    gta = VCFGenotypeAlleles.parse(r[sample], r.ref)
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
        return self.N() == self.gt_str.replace(".", "N").replace("|", "/")

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

    def alt_ploidy(self, allele=False):
        if allele:
            return len([x for x in self.allele_gt.split("/") if x == allele])
        else:
            return len([x for x in self.allele_gt.split("/") if x != self.ref])

    def is_ref(self):
        return len(set(self.alts) - set([self.ref])) == 0

    @staticmethod
    def parse(gt_with_allele_str, ref):
        alts = [x for x in gt_with_allele_str.split("/") if x != ref]
        alleles = [ref] + alts
        gt_str = "/".join([str(alleles.index(x)) for x in gt_with_allele_str.split("/")])
        gt = VCFGenotypeAlleles(gt_str, ref, alts)
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
    def create_gvcf(vcfs_folder, output_gvcf, ref_path, temp="/tmp/CombineVariants.vcf"):
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
        out_folder = os.path.dirname(temp)
        out_file = os.path.basename(temp) + ".bk"

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
        with open(temp) as h:
            if hasattr(output_gvcf, "write"):
                hw = output_gvcf
            else:
                hw = open(output_gvcf, "w")
            try:
                for l in h:
                    if l.startswith("#CHROM"):
                        vec = l.split("\t")
                        l = "\t".join(vec[:9] + [x.split(".variant")[0] for x in vec[9:]]) + "\n"

                    hw.write(l)
            finally:
                hw.close()

    @staticmethod
    def combineGVCFs(vcfs_folder, output_gvcf, ref_path, tmp="/tmp/combineGVCFs.vcf"):
        """

        :param vcfs_path_list: list of paths of the vcf files
        :param gvcf_path: gvcf to be created
        :param ref_path: fasta from the reference genome
        :return:
        """

        assert os.path.exists(ref_path), f'{ref_path} does no exists'
        assert os.path.exists(vcfs_folder), f'{vcfs_folder} does no exists'
        vcfs_folder = os.path.abspath(vcfs_folder)
        if not hasattr(output_gvcf, "write"):
            assert os.path.exists(os.path.dirname(output_gvcf)), f'{os.path.dirname(output_gvcf)} does no exists'

        vcf_files = []
        for x in glob(vcfs_folder + "/*vcf*"):
            if x.endswith(".vcf") or x.endswith(".vcf.gz") or x.endswith(".gvcf") or x.endswith(".gvcf.gz"):
                vcf_files.append(x)

        if not vcf_files:
            raise FileNotFoundError(f'no .vcf or .vcf.gz files where found at {vcfs_folder}')

        vcfs = " ".join([f"--variant {x}" for x in vcf_files])

        cmd = f"""
        gatk CombineGVCFs -R {ref_path} {vcfs} -O {tmp}
        """

        execute(cmd)

        with open(tmp) as h:
            if hasattr(output_gvcf, "write"):
                hw = output_gvcf
            else:
                hw = open(output_gvcf, "w")
            try:
                for l in h:
                    if l.startswith("#CHROM"):
                        vec = l.split("\t")
                        l = "\t".join(vec[:9] + [x.split(".variant")[0] for x in vec[9:]])
                    hw.write(l)
            finally:
                hw.close()

    @staticmethod
    def phylo(vcf, output):
        cmd = f"""bcftools filter -i 'alt=\"*\"' {vcf}  | bcftools norm -m -any | \
         bcftools filter -e 'alt=\"*\"'  | bcftools filter -i 'FORMAT/AD[*:1]>15' | \
         sed  's|0/1:|1/1|'  | sed  's|0\|1:|1/1|'  > /tmp/spaning_del.vcf"""
        execute(cmd)
        cmd = f"bcftools filter -e 'alt=\"*\"' {vcf}  > /tmp/no_spanning.vcf"
        execute(cmd)
        cmd = f"bcftools view /tmp/spaning_del.vcf | grep -v '^#' >> /tmp/no_spanning.vcf"
        execute(cmd)
        cmd = f"bcftools sort /tmp/no_spanning.vcf  > {output}"
        execute(cmd)

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
                    print(gt.N(), gt.isN(), gt)
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

        df = df_variants
        idx = df.pos != df.pos
        for i in range(len(samples)):
            s1 = df[samples[i]]
            s2 = df[samples[(i + 1) % len(samples)]]
            idx2 = s1 != s2
            idx = idx | idx2

        df = df[idx]

        return df[columns + samples]

    @staticmethod
    def aln(h, output, refseq=None, included_samples=None):

        # if hasattr(vcf_file, "read"):
        #     h = vcf_file
        # else:
        #     h = open(vcf_file)
        try:
            base_idx = 0
            for line in h:
                if line.startswith("#"):
                    if line.startswith("#CHROM"):
                        samples = [x.strip() for x in line.split()[9:]]

                        seqmap = {s: "" for s in samples}
                    continue
                break

            for line in tqdm(h, file=sys.stderr):
                _, pos, _, ref, alts = line.split()[:5]
                pos = int(pos)
                alts = [ref] + alts.split(",")
                gts = [x[0] for x in line.split()[9:]]
                gts = ["N" if x[0] == "." else alts[int(x[0])] for x in gts]
                pos_size = max([len(x) for x in alts])
                for i, s in enumerate(samples):
                    if not included_samples or s in included_samples:
                        subseq = refseq[base_idx:pos] + gts[i].ljust(pos_size, "-")
                        seqmap[s] += subseq

                sizes = {}
                for s in samples:
                    if not included_samples or s in included_samples:
                        sizes[s] = len(seqmap[s])
                assert len(set(sizes.values())) == 1, [base_idx, set(sizes.values()),
                                                       json.dumps({k: [x[0] for x in v] for k, v in
                                                                   groupby(sorted(sizes.items(), key=lambda x: x[1]),
                                                                           lambda x: x[1])})]

                base_idx = pos + len(ref)
        finally:
            h.close()

        for s in samples:
            if not samples or s in included_samples:
                seqmap[s] += refseq[base_idx:]

        if hasattr(output, "write"):
            h = output
        else:
            h = open(output, "w")
        try:
            for k, v in seqmap.items():
                bpio.write(SeqRecord(id=k, name="", description="", seq=Seq(v)), h, "fasta")
        finally:
            h.close()


if __name__ == '__main__':
    import argparse
    from SNDG import init_log
    import fileinput

    init_log()

    parser = argparse.ArgumentParser(description='Utils for vcfs with multiple samples')

    subparsers = parser.add_subparsers(help='commands', description='valid subcommands', dest='command', required=True)
    cmd = subparsers.add_parser('combine', help='annotates a list of residues')
    cmd.add_argument('-i', '--vcfs_dir', required=True, help='directory with vcfs to combine')
    cmd.add_argument('-ref', '--reference', required=True, help='fasta reference')

    cmd = subparsers.add_parser('phylo', help='prepare vcf for phylogeny')
    cmd.add_argument('-i', '--vcf', required=True, help='joint vcf file')
    cmd.add_argument('-ref', '--reference', required=True, help='fasta reference')

    cmd = subparsers.add_parser('aln', help='fasta aln')
    cmd.add_argument('-i', '--vcf', default="-", help='joint vcf file')
    cmd.add_argument('-ref', '--reference', required=True, help='fasta reference')
    cmd.add_argument('--include', nargs='*', default=[])

    args = parser.parse_args()

    if args.command == "combine":
        VariantSetUtils.combineGVCFs(args.vcfs_dir, sys.stdout, args.reference)
    if args.command == "phylo":
        VariantSetUtils.phylo(args.vcf, sys.stdout)
    if args.command == "aln":
        if args.include and os.path.exists(args.include[0]):
            sys.stderr.write(f"importing list from {args.include[0]}\n")
            samples = [x.strip() for x in open(args.include[0]).readlines() if x.strip()]
        else:
            samples = args.include
        if samples:
            sys.stderr.write(f'filtering {len(samples)} samples:{",".join(samples)}\n')

        VariantSetUtils.aln(fileinput.input(args.vcf), sys.stdout,
                            refseq=str(bpio.read(args.reference, "fasta").seq),
                            included_samples=samples)

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
    #
    # def default_value_fn(variant, sample, alt, assigned_values):
    #     if sample.sample == "3867INF":
    #         s = [x for x in variant.samples if x.sample == "3867NI"][0]
    #         if s.called and validate_variant_fn(s.data):
    #             alt = str(variant.ALT[int(s.data.GT) - 1])
    #             return alt
    #         else:
    #             s = [x for x in variant.samples if x.sample == "3867NE"][0]
    #             if s.called and validate_variant_fn(s.data):
    #                 alt = str(variant.ALT[int(s.data.GT) - 1])
    #                 return alt
    #             return variant.REF
    #     elif sample.sample == "3867NE":
    #         s = [x for x in variant.samples if x.sample == "3867INF"][0]
    #         if s.called and validate_variant_fn(s.data):
    #             alt = str(variant.ALT[int(s.data.GT) - 1])
    #             return alt
    #         else:
    #             s = [x for x in variant.samples if x.sample == "3867NI"][0]
    #             if s.called and validate_variant_fn(s.data):
    #                 alt = str(variant.ALT[int(s.data.GT) - 1])
    #                 return alt
    #             return variant.REF
    #     elif sample.sample == "3867NI":
    #         s = [x for x in variant.samples if x.sample == "3867INF"][0]
    #         if s.called and validate_variant_fn(s.data):
    #             return str(variant.ALT[int(s.data.GT) - 1])
    #         else:
    #             s = [x for x in variant.samples if x.sample == "3867INF"][0]
    #             if s.called and validate_variant_fn(s.data):
    #                 return str(variant.ALT[int(s.data.GT) - 1])
    #             return variant.REF
    #     return alt

    # vs = VariantSet("/mnt/data2/data/projects/23staphylo_old/external/sordelli/variant_call/cohort.g.vcf")
    # # vs.default_not_called_fn = default_not_called_fn
    # vs.default_value_fn = default_value_fn
    #
    # df = vs.build_table()
    # with open("/tmp/imputed.vcf", "w") as h:
    #     VariantSet.complete_vcf_with_bams("/home/eze/projects/ST5/processed/groups/arg_chi/all.vcf",
    #                                       "/home/eze/projects/ST5/processed/strains/*/alignment/mapped_reads.bam",
    #                                       stdout=h)
