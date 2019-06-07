"""

"""
import os
import subprocess as sp

import pandas as pd
from tqdm import tqdm
from collections import defaultdict

from SNDG import execute
from SNDG.Comparative.VcfSnpeffIO import VcfSnpeffIO
import pysam
import Bio.SeqIO as bpio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def validate_variant_fn(sample_data,min_depth=30):
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


class VariantSet():
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
    def create_gvcf(vcfs_path_list, gvcf_path, ref_path):
        """

        :param vcfs_path_list: list of paths of the vcf files
        :param gvcf_path: gvcf to be created
        :param ref_path: fasta from the reference genome
        :return:
        """
        cmd_template = """
        java -jar $GATK -T CombineVariants    -R {ref} {vcfs}   -o {out}    -genotypeMergeOptions UNIQUIFY
        """
        vcfs = " ".join(["--variant " + x for x in vcfs_path_list])
        cmd = cmd_template.format(vcfs=vcfs, out=gvcf_path, ref=ref_path)
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

    def build_table(self, keep_all=False,min_depth=30):
        class Empty:
            pass

        rows = []
        for variant, effects in tqdm(self.gvcf, total=self.total_variants):
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
                "pos": variant.POS - 1,
                "impact": effect.impact,
                "type": "&".join(effect.annotation)}
            alternatives = [variant.REF] if keep_all else []
            for sample in variant.samples:

                sample_name = sample.sample.split(".variant")[0]

                if sample.called:
                    assert sample.data.GT != "."
                    alt = str(variant.ALT[int(sample.data.GT) - 1])
                    if alt == "*":
                        alt = ""

                    if hasattr(sample.data, "AD") and isinstance(sample.data.AD, list):
                        vresult[sample_name + "_ada"] = sample.data.AD[1]
                        vresult[sample_name + "_adr"] = sample.data.AD[0]
                    elif hasattr(sample.data, "AD"):
                        vresult[sample_name + "_ada"] = sample.data.AD
                        vresult[sample_name + "_adr"] = 0
                    else:
                        vresult[sample_name + "_ada"] = 0
                        vresult[sample_name + "_adr"] = 0

                    if self.validate_variant_fn(sample.data,min_depth):
                        vresult[sample_name] = alt
                    else:
                        vresult[sample_name] = self.default_value_fn(variant, sample, alt, dict(vresult))
                    if (alt == vresult[sample_name]) and effect.aa_pos:
                        vresult["aa_pos"] = effect.aa_pos
                        vresult["aa_ref"] = effect.aa_ref
                        vresult["aa_alt"] = effect.aa_alt
                elif sample.data.GT == ".":
                    vresult[sample_name] = self.default_value_fn(variant, sample, variant.REF, dict(vresult))
                else:
                    vresult[sample_name] = variant.REF
                    if self.bam:
                        for pileupcolumn in self.bam[sample_name].pileup(variant.CHROM, variant.POS - 1, variant.POS,
                                                                         truncate=True):  # min_base_quality = 0, stepper = "nofilter"
                            pileup = []
                            indel_pileup = []
                            for pileupread in pileupcolumn.pileups:
                                if not pileupread.is_del and not pileupread.is_refskip:
                                    # query position is None if is_del or is_refskip is set.
                                    pileup.append(pileupread.alignment.query_sequence[pileupread.query_position])
                                    indel_pileup.append(pileupread.indel)
                        if 'ins' in str(effect.hgvs_c) or 'dup' in str(effect.hgvs_c) or 'del' in str(effect.hgvs_c):
                            cov_alt = len([x for x in indel_pileup if x != 0])
                            cov_ref = len(indel_pileup) - cov_alt
                        else:
                            cov_alt = len([x for x in pileup if x != variant.REF])
                            cov_ref = len(pileup) - cov_alt
                        vresult[sample_name + "_ada"] = cov_alt
                        vresult[sample_name + "_adr"] = cov_ref
                        if self.default_not_called_fn(cov_ref, cov_alt):
                            vresult[sample_name] = self.default_value_fn(variant, sample, alt, dict(vresult))
                        if ('ins' in str(effect.hgvs_c) or 'dup' in str(effect.hgvs_c)) and len(
                                max(variant.ALT, key=len)) >= 25:
                            # hay problemas en el bam con los ins muy largos
                            vresult[sample_name] = self.default_value_fn(variant, sample, alt, dict(vresult))

                alternatives.append(vresult[sample_name])

            if len(set(alternatives)) > 1:
                rows.append(vresult)
                # df_result = df_result.append(vresult, ignore_index=True)
        if self.bam:
            for sample in self.bam:
                self.bam[sample].close()

        df = pd.DataFrame(rows)

        return df

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

    def aln(self,aln,txt="/tmp/samples.txt",ref=None):

        ref = str(bpio.read(ref,"fasta").seq)
        if not os.path.exists(txt):
            df = self.build_table()
            df.to_csv(txt)
        else:
            df = pd.read_csv(txt)
            df = df.fillna("")
        samples = [c for c in df.columns if c not in ["chrom","pos","ref",'Unnamed: 0',"gene","impact","type",]
                   and  ("_ada" not in c) and ("_adr" not in c)]

        seqmap = {s: "" for s in samples}
        total = len(df)
        base_idx = 0
        for _, r in tqdm(df.sort_values("pos").iterrows(), total=total):
            pos_size = max([len(r[x]) for x in samples ])
            for s in samples:
                subseq = ref[base_idx:r.pos] +  r[s].ljust(pos_size, "-")
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


    vs = VariantSet("/mnt/data2/data/projects/23staphylo_old/external/sordelli/variant_call/cohort.g.vcf")
    # vs.default_not_called_fn = default_not_called_fn
    vs.default_value_fn = default_value_fn

    df = vs.build_table()
