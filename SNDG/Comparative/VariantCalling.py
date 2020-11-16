import multiprocessing
import os

import pandas as pd
from tqdm import tqdm

from SNDG import execute as e, docker_wrap_command
from SNDG import mkdir
from SNDG.Comparative.CoverageAnalysis import CoverageAnalysis


class VariantCalling:

    @staticmethod
    def genotype_call(reference, vcf, output_file="./combined.vcf", ploidy=2):
        wd = os.path.dirname(os.path.abspath(output_file)) + "/"
        reference = os.path.abspath(reference)
        vcf = os.path.abspath(vcf)

        mkdir(wd)

        assert os.path.exists(wd), f'{wd} could not be created'
        assert os.path.exists(reference), f'{reference} does not exist'
        assert os.path.exists(vcf), f'{vcf} does not exist'

        e(f"""gatk GenotypeGVCFs \
        -R "{reference}" -ploidy {ploidy} \
        -V "{vcf}" \
        -O "{output_file}" 
        """)
        return

    @staticmethod
    def variant_call(wd, reference, alignment, ploidy=2):
        wd = os.path.abspath(wd) + "/"
        reference = os.path.abspath(reference)
        alignment = os.path.abspath(alignment)

        mkdir(wd)

        assert os.path.exists(wd), f'{wd} could not be created'
        assert os.path.exists(reference), f'{reference} does not exist'
        assert os.path.exists(alignment), f'{alignment} does not exist'

        e(f"""gatk HaplotypeCaller -ERC GVCF \
         -R "{reference}" -ploidy {ploidy} \
         -I "{alignment}" --output-mode EMIT_ALL_CONFIDENT_SITES \
         -O "{wd}raw.g.vcf.gz"
        """)

        e(f"""gatk GenotypeGVCFs \
        -R "{reference}" -ploidy {ploidy} \
        -V "{wd}raw.g.vcf.gz" \
        -O "{wd}output.vcf.gz" 
        """)

        # # Call variants in the sequence data
        # e(
        #     "java -jar $GATK -T HaplotypeCaller -R {record_name} -I {alignment} -gt_mode DISCOVERY -ploidy 1 -stand_call_conf 30 -o raw_variants.vcf"
        #     , work_dir, record_name=record, alignment=alignment)
        # # Apply hard filters to a call set
        # e("java -jar $GATK -T SelectVariants -R {record_name} -V raw_variants.vcf  -selectType SNP -o raw_snps.vcf"
        #   , work_dir, record_name=record)
        # e(
        #     "java -jar $GATK -T VariantFiltration -R {record_name} -V raw_snps.vcf  -filter \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" --filterName \"my_snp_filter\" -o filtered_snps.vcf",
        #     work_dir, record_name=record)
        # e("java -jar $GATK -T SelectVariants -R {record_name} -V raw_variants.vcf  -selectType INDEL -o raw_indels.vcf",
        #   work_dir, record_name=record)
        # e(
        #     "java -jar $GATK -T VariantFiltration -R {record_name} -V raw_indels.vcf  -filter \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\" --filterName \"my_indel_filter\" -o filtered_indels.vcf",
        #     work_dir, record_name=record)
        # e(
        #     "java -jar $GATK -T CombineVariants --assumeIdenticalSamples -R {record_name} -V filtered_snps.vcf -V filtered_indels.vcf -genotypeMergeOptions UNIQUIFY -o concatenated.vcf",
        #     work_dir, record_name=record)
        # # Removes column from vcf header
        # e("sed \'/^#[^#]/ {{s/\\t%s\\.variant2//}}\' concatenated.vcf > %s.vcf" % (strain, "final.vcf"), work_dir)
        # return strain + ".vcf"

    # @staticmethod
    # def variant_call(in_bam, ref_fasta, snpeffdb, work_dir="./", cpus=1):
    #     mkdir(work_dir)
    #     variants_file = work_dir + "/variants.vcf"
    #     variants_ann_file = work_dir + "/variants.ann.vcf"
    #
    #     if not os.path.exists(variants_file):
    #         e(
    #             "gatk -T HaplotypeCaller -ploidy 1 -R {ref} -I {aln} --num_cpu_threads_per_data_thread {cpus} --genotyping_mode DISCOVERY -stand_call_conf 30 -o {vcf}",
    #             ref=ref_fasta, aln=in_bam, vcf=variants_file, cpus=cpus)
    #     if not os.path.exists(variants_ann_file):
    #         e("snpEff {database} {input} > {output}",
    #                 database=snpeffdb, input=variants_file, output=variants_ann_file)

    @staticmethod
    def create_snpeff_db():
        """
        cat NC_002505.1.gbk NC_002506.1.gbk > genes.gbk
        vibrio.genome : Vibrio Cholerae
	        vibrio.chromosomes : NC_002505.1, NC_002506.1
	        vibrio.NC_002505.1.codonTable : Bacterial_and_Plant_Plastid
	    java -jar snpEff.jar build -genbank -v CP000730
        :return:
        """
        pass

    @staticmethod
    def vcf_stats(in_vcf):
        pass


if __name__ == '__main__':
    import argparse
    from SNDG import init_log

    init_log()

    parser = argparse.ArgumentParser(description='Mapping to variant calls pipeline.')

    subparsers = parser.add_subparsers(help='commands', description='valid subcommands', dest='command', required=True)
    cmd = subparsers.add_parser('vc', help='variant calling pipeline')
    cmd.add_argument('-ref', '--reference', required=True, help='fasta file. Can be gziped')
    cmd.add_argument('-b', '--bam', required=True)
    cmd.add_argument('-w', '--workdir', default="./")
    cmd.add_argument('-p', '--ploidy', default=2, type=int)

    cmd = subparsers.add_parser('gt', help='genotyping for multi sample')
    cmd.add_argument('-ref', '--reference', required=True, help='fasta file. Can be gziped')
    cmd.add_argument('--vcf', required=True)
    cmd.add_argument('-o', '--ouput_vcf', default="./output.vcf")
    cmd.add_argument('-p', '--ploidy', default=2, type=int)

    cmd = subparsers.add_parser('aa', help='annotate vcf')
    cmd.add_argument('-w', '--workdir', default="./")
    cmd.add_argument('--cpus', type=int, default=multiprocessing.cpu_count())

    # parser.add_argument('-s', '--steps', nargs='*', default=["binding", "pocket"], choices=["binding", "pocket"])

    # parser.add_argument('--useSingletons', action = 'store_true', dest = 'singletons')

    args = parser.parse_args()
    if args.command == 'vc':
        VariantCalling.variant_call(args.workdir, args.reference, args.bam,ploidy=args.ploidy)
    if args.command == 'gt':
        VariantCalling.genotype_call( args.reference, args.vcf,args.ouput_vcf,ploidy=args.ploidy)
    if args.command == 'ann':
        raise NotImplemented()
