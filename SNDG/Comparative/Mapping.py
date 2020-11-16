import multiprocessing
import os

import pandas as pd
from tqdm import tqdm

from SNDG import execute as e, docker_wrap_command
from SNDG import mkdir
from SNDG.Comparative.CoverageAnalysis import CoverageAnalysis


class Mapping:

    @staticmethod
    def init_ref(reference_path):
        last = reference_path.split(".")[-1]
        if last == "gz":
            last = reference_path.split(".")[-2] + ".gz"
        assert last in ["fna", "fasta", "fa", "fna.gz", "fasta.gz", "fa.gz"], f'unknown extension for {reference_path}'
        dict_path = reference_path.replace("." + last, ".dict")
        e(f"bwa index -a is  {reference_path}")
        e(f"samtools dict {reference_path} > {dict_path}")
        if last.endswith(".gz"):
            e(
                f"zcat {reference_path}| bgzip > {reference_path}.tmp; cp '{reference_path}.tmp' '{reference_path}' && rm '{reference_path}.tmp'")

        e(f"samtools faidx {reference_path}")
        # e("bowtie2-build {record_name} {record_name}".format(record_name=filename))
        if reference_path.endswith(".gz"):
            e(
                f"zcat {reference_path} | makeblastdb -dbtype nucl -title {reference_path} -input_type fasta -out {reference_path}  -in -")
        else:
            e(f"makeblastdb -dbtype nucl -in {reference_path} ")

    @staticmethod
    def clean_reads(work_dir, read1, read2, trim_left=20,
                    trim_qual_right=25, trim_qual_window=25, min_len=35, window_size=5, cpu=1,
                    clip_string=""):
        """

        :param work_dir:
        :param read1:
        :param read2:
        :param min_qual_mean:
        :param trim_left:
        :param trim_qual_right:
        :param trim_qual_window:
        :param min_len:
        :param window_size:
        :param cpu:
        :param clip_string: example -> ILLUMINACLIP:TruSeq3-SE:2:30:10
        :return:
        """
        work_dir = os.path.abspath(work_dir) + "/"
        read1 = os.path.abspath(read1)
        read2 = os.path.abspath(read2)

        # Quality control
        # "prinseq-lite.pl -fastq {read1_full} -fastq2 {read2_full} -min_qual_mean {min_qual_mean}" +
        # " -trim_left {trim_left}  -trim_qual_right {trim_qual_right} -trim_qual_window {trim_qual_window}" +
        # " -min_len {min_len} -out_good trimmed",
        e(
            'java -jar $TRIMMOMATIC PE -threads {cpu} "{read1_full}" "{read2_full}" ' +
            ' {pout1} {upout1} {pout2} {upout2} ' +
            ' HEADCROP:{trim_left}  TRAILING:{trim_qual_right} SLIDINGWINDOW:{window_size}:{trim_qual_window} ' +
            ' {clip_string}  MINLEN:{min_len} ',
            work_dir, read1_full=read1, read2_full=read2, trim_left=trim_left, trim_qual_right=trim_qual_right,
            trim_qual_window=trim_qual_window, min_len=min_len, cpu=cpu,
            pout1="trimmed_1.fastq", pout2="trimmed_2.fastq", upout1="trimmed_1_singletons.fastq",
            upout2="trimmed_2_singletons.fastq", window_size=window_size, clip_string=clip_string
        )

        e("fastqc trimmed_1.fastq", work_dir)
        e("fastqc trimmed_2.fastq", work_dir)

        if os.path.exists(work_dir + "trimmed_1_singletons.fastq"):
            e("cat trimmed_1_singletons.fastq >> trimmed_s.fastq", work_dir)
            os.remove(work_dir + "trimmed_1_singletons.fastq")

        if os.path.exists(work_dir + "trimmed_2_singletons.fastq"):
            e("cat trimmed_2_singletons.fastq >> trimmed_s.fastq", work_dir)
            os.remove(work_dir + "trimmed_2_singletons.fastq")

    @staticmethod
    def alignment(wd, ref, trimmed_1="trimmed_1.fastq",
                  trimmed_2="trimmed_2.fastq", cpus=multiprocessing.cpu_count(),
                  strain="sample1", species=None,force=False,read_group="group1"):
        if not species:
            species = strain

        mkdir(wd)

        wd = os.path.abspath(wd) + "/"
        ref = os.path.abspath(ref)

        assert os.path.exists(wd),f'{wd} could not be created'
        assert os.path.exists(ref),f'{ref} does not exist'
        assert os.path.exists(trimmed_1),f'{trimmed_1} does not exist'
        assert os.path.exists(trimmed_2),f'{trimmed_2} does not exist'

        # Generate a SAM file containing aligned reads
        if force or not os.path.exists(f"{wd}mapped_reads_raw.bam"):
            tab = "\\t"
            e(
                f"bwa mem -t {cpus} -M -R \'@RG{tab}ID:{read_group}{tab}SM:{strain}{tab}PL:illumina{tab}LB:{species}\' {ref} {trimmed_1} {trimmed_2} > {wd}aligned_reads.sam")
        assert os.path.getsize(f"{wd}aligned_reads.sam") > 10, f"{wd}aligned_reads.sam cant be empty"
        # Filter mapped reads and convert to BAM
        if force or (not os.path.exists(f"{wd}dedup.bam ") and not os.path.exists(f"{wd}mapped_reads_raw.bam")):
            e(f"samtools view -@ {cpus} -F 4 -S -b -h {wd}aligned_reads.sam | samtools sort - > {wd}mapped_reads_raw.bam")
            e(f"samtools view -@ {cpus} -f 4 -S -b -h {wd}aligned_reads.sam > {wd}unmapped_reads.bam")
            e(f"bedtools bamtofastq -i unmapped_reads.bam -fq {wd}unmapped_1.fastq -fq2 {wd}unmapped_2.fastq")
        if os.path.exists(f"{wd}unmapped_reads.bam"):
            os.remove(f"{wd}unmapped_reads.bam")

        if os.path.exists(f"{wd}aligned_reads.sam"):
            os.remove(f"{wd}aligned_reads.sam")

        # Sort and mark duplicates
        e(f"gatk MarkDuplicates -INPUT {wd}mapped_reads_raw.bam -OUTPUT {wd}dedup.bam -METRICS_FILE {wd}metrics.txt")
        assert os.path.getsize(f"{wd}dedup.bam") > 10, f"{wd}dedup.bam cant be empty"

        os.remove(f"{wd}mapped_reads_raw.bam")
        e(f'samtools sort {wd}dedup.bam > {wd}mapped_reads.bam')
        os.remove(f"{wd}dedup.bam")
        e(f'samtools index {wd}mapped_reads.bam')
        e(
            f"gatk CollectInsertSizeMetrics --I {wd}mapped_reads.bam --O {wd}insert_size_metrics.txt --H {wd}insert_size_histogram.pdf --M 0.5")

        return f'{wd}mapped_reads.bam'



    @staticmethod
    def realign(bam_file, ref_fasta):
        out_bwa_bam = "sorted_" + bam_file
        e("samtools sort -o %s %s" % (out_bwa_bam, bam_file))
        out_bwa_final_bam = "realigned2_" + bam_file
        out_bwa_intervals = bam_file + ".intervals"
        out_bwa_intervals2 = bam_file + ".intervals"
        bwa_realigned = "realigned_" + bam_file
        duplicates = "duplicates_" + bam_file
        bwa_iter1 = "iter1_" + bam_file
        if not os.path.exists(out_bwa_final_bam):
            e("samtools index %s" % out_bwa_bam)
            e("gatk -T RealignerTargetCreator -R {ref} -I {input} -o {out}",
              ref=ref_fasta, input=out_bwa_bam, out=out_bwa_intervals)
            e("gatk -T IndelRealigner -R {ref} -I {input} -targetIntervals {intervals} -o {output}",
              ref=ref_fasta, input=out_bwa_bam, intervals=out_bwa_intervals, output=bwa_realigned)
            # Aca se recomienda correr el BaseRecalibrator de GATK pero no se tiene un vcf con variantes comunes
            e("picard MarkDuplicates I={input}   REMOVE_DUPLICATES=true O={output} M={duplicates}",
              input=bwa_realigned, output=bwa_iter1, duplicates=duplicates)
            e("samtools index {input}", input=bwa_iter1)
            e("gatk -T RealignerTargetCreator -R {ref} -I {input} -o {intervals}",
              ref=ref_fasta, input=bwa_iter1, intervals=out_bwa_intervals2)
            e("gatk -T IndelRealigner -R {ref} -I {input} -targetIntervals {intervals} -o {output}",
              ref=ref_fasta, input=bwa_iter1, intervals=out_bwa_intervals2, output=out_bwa_final_bam)
            e("samtools index %s" % out_bwa_final_bam)

        for x in [bam_file, out_bwa_intervals, out_bwa_intervals2, bwa_realigned, bwa_iter1, bwa_iter1 + ".bai"]:
            if os.path.exists(x):
                os.remove(x)

        return out_bwa_final_bam



if __name__ == '__main__':
    import argparse
    from SNDG import init_log

    init_log()

    parser = argparse.ArgumentParser(description='Mapping to variant calls pipeline.')

    subparsers = parser.add_subparsers(help='commands', description='valid subcommands', dest='command', required=True)
    cmd = subparsers.add_parser('ref_init', help='annotates a list of residues')
    cmd.add_argument('-i', '--reference', required=True, help='fasta file. Can be gziped')

    cmd = subparsers.add_parser('aln', help='runs bam creation pipeline')
    cmd.add_argument('-w', '--workdir', default="./")
    cmd.add_argument('-ref', '--reference', required=True)
    cmd.add_argument('-s', '--sample', default="sample")
    cmd.add_argument('-x', '--species', default="-")
    cmd.add_argument('-r1', '--read1', required=True)
    cmd.add_argument('-r2', '--read2', required=True)
    cmd.add_argument( '--force', action="store_true")
    cmd.add_argument('--cpus', type=int, default=multiprocessing.cpu_count())

    # parser.add_argument('-s', '--steps', nargs='*', default=["binding", "pocket"], choices=["binding", "pocket"])

    # parser.add_argument('--useSingletons', action = 'store_true', dest = 'singletons')

    args = parser.parse_args()
    if args.command == 'ref_init':
        Mapping.init_ref(args.reference)
    if args.command == 'aln':
        Mapping.alignment(args.workdir, trimmed_1=args.read1, ref=args.reference,
                          trimmed_2=args.read1, cpus=args.cpus,
                          strain=args.sample, species=args.species,force=args.force)


