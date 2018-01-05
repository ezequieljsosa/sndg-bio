import os
from SNDG import execute

"""
requires:
samtools bcftools bowtie bedtools and bwa binaries in PATH
ncbi-blast+ installed
picard and gatk link in path:
example:
/usr/local/bin/gatk
#!/bin/bash
java -jar /opt/GATK/GenomeAnalysisTK.jar $@



"""

def init_ref(path, record):
    cwd = os.getcwd()
    try:
        os.chdir(path)
        execute("bwa index -a is  {record_name}.fasta".format(record_name=record))
        execute("picard CreateSequenceDictionary R={record_name}.fasta O={record_name}.dict".format(record_name=record))
        execute("samtools faidx {record_name}.fasta".format(record_name=record))
        execute("bowtie2-build {record_name}.fasta {record_name}".format(record_name=record))
        execute("makeblastdb -dbtype nucl -in {record_name}.fasta ".format(record_name=record))
    finally:
        os.chdir(cwd)

def realign(bam_file, ref_fasta):
    out_bwa_bam = "sorted_" + bam_file
    execute("samtools sort -o %s %s" % (out_bwa_bam, bam_file))
    out_bwa_final_bam = "realigned2_" + bam_file
    out_bwa_intervals = bam_file + ".intervals"
    out_bwa_intervals2 = bam_file + ".intervals"
    bwa_realigned = "realigned_" + bam_file
    duplicates = "duplicates_" + bam_file
    bwa_iter1 = "iter1_" + bam_file
    if not os.path.exists(out_bwa_final_bam):
        execute("samtools index %s" % out_bwa_bam)
        execute("gatk -T RealignerTargetCreator -R {ref} -I {input} -o {out}",
                ref=ref_fasta, input=out_bwa_bam, out=out_bwa_intervals)
        execute("gatk -T IndelRealigner -R {ref} -I {input} -targetIntervals {intervals} -o {output}",
                ref=ref_fasta, input=out_bwa_bam, intervals=out_bwa_intervals, output=bwa_realigned)
        # Aca se recomienda correr el BaseRecalibrator de GATK pero no se tiene un vcf con variantes comunes
        execute("picard MarkDuplicates I={input}   REMOVE_DUPLICATES=true O={output} M={duplicates}",
                input=bwa_realigned, output=bwa_iter1, duplicates=duplicates)
        execute("samtools index {input}", input=bwa_iter1)
        execute("gatk -T RealignerTargetCreator -R {ref} -I {input} -o {intervals}",
                ref=ref_fasta, input=bwa_iter1, intervals=out_bwa_intervals2)
        execute("gatk -T IndelRealigner -R {ref} -I {input} -targetIntervals {intervals} -o {output}",
                ref=ref_fasta, input=bwa_iter1, intervals=out_bwa_intervals2, output=out_bwa_final_bam)
        execute("samtools index %s" % out_bwa_final_bam)

    for x in [bam_file, out_bwa_intervals, out_bwa_intervals2, bwa_realigned, bwa_iter1, bwa_iter1 + ".bai"]:
        if os.path.exists(x):
            os.remove(x)

    return out_bwa_final_bam


def process_ref(ref_fasta, strains, read_template, work_dir):
    """

    read_template example: "/data/projects/Staphylococcus/samples/trimmed/%s.%s.fastq" -> path + (strain,R1|R2|S)
    """
    print "oh yeah!2"
    out_bwa_bam = "final_bwa.bam"
    out_bwa_bam_idx = out_bwa_bam + ".bai"

    for cepa in strains:

        if not os.path.exists(work_dir + cepa):
            os.makedirs(work_dir + cepa)
        os.chdir(work_dir + cepa)

        out_bwa_pe = "bwa_pe.sam"
        out_bwa_pe_bam = "bwa_pe.bam"
        out_unmapped_pe_bam = "unmapped.bam"

        if not os.path.exists(out_bwa_bam_idx) and not os.path.exists(out_bwa_pe_bam):
            execute('bwa mem -R "@RG\\tID:illumina\\tSM:{ncepa}\\tLB:{ncepa}"  {ref_fasta} {pe1} {pe2}  >  ' + out_bwa_pe,
                    ref_fasta=ref_fasta, ncepa=cepa,
                    pe1=read_template % (cepa, "R1"),
                    pe2=read_template % (cepa, "R2"))
            execute("samtools view -F 4 -Sbh %s > %s" % (out_bwa_pe, out_bwa_pe_bam))
            execute("samtools view -f 4 -Sbh %s > %s" % (out_bwa_pe, out_unmapped_pe_bam))
            unmapped_pair_1 = "unmapped_pair_1.fastq"
            unmapped_pair_2 = "unmapped_pair_2.fastq"
            execute("bedtools bamtofastq -i {ubam} -fq {upair}   -fq2 {upair2}",
                    upair2=unmapped_pair_2, upair=unmapped_pair_1, ubam=out_unmapped_pe_bam)
            out_bwa_pe_bam = realign(out_bwa_pe_bam, ref_fasta)

        for x in [out_bwa_pe, out_unmapped_pe_bam]:
            if os.path.exists(x):
                os.remove(x)

        out_bwa_se = "bwa_se.sam"
        out_bwa_se_bam = "bwa_se.bam"
        out_unmapped_se_bam = "unmapped_se.bam"
        if not os.path.exists(out_bwa_bam_idx) and not os.path.exists(out_bwa_se_bam):
            execute('bwa mem -R "@RG\\tID:illumina\\tSM:{ncepa}\\tLB:{ncepa}"  {ref_fasta} {s1}   >  ' + out_bwa_se,
                    ref_fasta=ref_fasta, ncepa=cepa,
                    s1=read_template % (cepa, "S"))
            execute("samtools view -F 4 -Sbh %s > %s" % (out_bwa_se, out_bwa_se_bam))
            execute("samtools view -f 4 -Sbh %s > %s" % (out_bwa_se, out_unmapped_se_bam))

            unmapped_single = "unmapped_single.fastq"
            execute("bedtools bamtofastq -i {ubam} -fq {upair}   ",
                    upair=unmapped_single, ubam=out_unmapped_se_bam)
            out_bwa_se_bam = realign(out_bwa_se_bam, ref_fasta)

        for x in [out_bwa_se, out_unmapped_se_bam]:
            if os.path.exists(x):
                os.remove(x)

        out_bwa_raw_bam = "bwa_raw.bam"
        out_bwa_fm_bam = "bwa_fm.bam"
        out_bwa_fm_sort_bam = "bwa_fm_sort.bam"

        if not os.path.exists(out_bwa_bam_idx):
            execute("samtools merge %s %s %s " % (out_bwa_raw_bam, out_bwa_pe_bam, out_bwa_se_bam))
            execute("samtools sort -n -o %s %s" % (out_bwa_fm_sort_bam, out_bwa_raw_bam))
            execute("samtools fixmate  %s %s" % (out_bwa_fm_sort_bam, out_bwa_fm_bam))
            execute("samtools sort -o %s %s" % (out_bwa_bam, out_bwa_fm_bam))
            execute("samtools index %s" % out_bwa_bam)

        for x in [out_bwa_raw_bam, out_bwa_fm_bam, out_bwa_fm_sort_bam,
                  out_bwa_pe_bam, out_bwa_se_bam, out_bwa_pe_bam + ".bai", out_bwa_se_bam + ".bai"]:
            if os.path.exists(x):
                os.remove(x)

        if not os.path.exists("flagstat.txt"):
            execute("samtools flagstat %s > %s" % (out_bwa_bam, "flagstat.txt"))

def variant_call(strain,  in_bam, ref_fasta, snpeffdb,work_dir="./"):
    if not os.path.exists(work_dir + strain):
        os.makedirs(work_dir + strain)
    os.chdir(work_dir + strain)
    variants_file = "variants.vcf"
    variants_ann_file = "variants.ann.vcf"



    if not os.path.exists(variants_file):
        execute("gatk -T HaplotypeCaller -ploidy 1 -R {ref} -I {aln} --num_cpu_threads_per_data_thread 4 --genotyping_mode DISCOVERY -stand_call_conf 30 -o {vcf}",
                ref=ref_fasta, aln=in_bam, vcf=variants_file)
    if not os.path.exists(variants_ann_file):
        execute("snpEff {database} {input} > {output}",
                database=snpeffdb, input=variants_file, output=variants_ann_file)