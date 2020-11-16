import os
import subprocess as sp
import shutil
from subprocess import CalledProcessError

def execute(cmd_unformated, wd="./", retcodes=[0], docker_mode=False, **kargs):
    cmd = "cd " + wd + ";" + cmd_unformated.format(**kargs)

    try:
        sp.check_output(cmd, shell=True, stderr=sp.STDOUT)
    except CalledProcessError as ex:
        if ex.returncode not in retcodes:
            raise

class MultiRefCoreAln:

    def __init__(self, base_refs, cpus=1):
        self.cpus = cpus
        self.base_refs = base_refs

    def core_reads(self, output, cepa,
                   in_fasta_r1, in_fasta_r2, in_fasta_s,
                   tmp_dir="/tmp/"):
        cepa_reads_core_r1_bk = output + cepa + "_R1.fastq.bk"
        cepa_reads_core_r2_bk = output + cepa + "_R2.fastq.bk"
        cepa_reads_core_s_bk = output + cepa + "_S.fastq.bk"

        cepa_reads_core_r1 = output + cepa + "_R1.fastq"
        cepa_reads_core_r2 = output + cepa + "_R2.fastq"
        cepa_reads_core_s = output + cepa + "_S.fastq"

        aln_pe_sam = f"{tmp_dir}/aln_re.sam"
        aln_pe_bam = f"{tmp_dir}/aln_re.bam"
        aln_s_sam = f"{tmp_dir}/aln_s.sam"
        aln_s_bam = f"{tmp_dir}/aln_s.bam"
        for ref in self.base_refs:
            ref_fasta = ref
            execute('bwa mem -t {cpus} -R "@RG\\tID:illumina\\tSM:{ncepa}\\tLB:{ncepa}"  {ref_fasta} {pe1} {pe2}  > ' +
                    aln_pe_sam, ref_fasta=ref_fasta, ncepa=cepa, pe1=in_fasta_r1, pe2=in_fasta_r2, cpus=self.cpus)
            execute("samtools view -@ %s  -F 4 -Sbh %s > %s" % (self.cpus, aln_pe_sam, aln_pe_bam))
            os.remove(aln_pe_sam)
            execute("bedtools bamtofastq -i {ubam} -fq {upair}   -fq2 {upair2}",
                    upair2=cepa_reads_core_r1, upair=cepa_reads_core_r2, ubam=aln_pe_bam)
            if os.path.getsize(cepa_reads_core_r1) < 1000:
                print("NOOOOOOOOO!!!!!!!!!!!!!!!!!!!: " + cepa + " -> " + ref)
            if in_fasta_s:
                execute(
                    'bwa mem -t {cpus} -R "@RG\\tID:illumina\\tSM:{ncepa}\\tLB:{ncepa}"  {ref_fasta} {s1}   >  ' + aln_s_sam,
                    ref_fasta=ref_fasta, ncepa=cepa, s1=in_fasta_s, cpus=self.cpus)
                execute("samtools view -@ %s -F 4 -Sbh %s > %s" % (self.cpus, aln_s_sam, aln_s_bam))
                os.remove(aln_s_sam)
                execute("bedtools bamtofastq -i {ubam} -fq {upair} ", upair=cepa_reads_core_s, ubam=aln_s_bam)
                in_fasta_s = cepa_reads_core_s
                shutil.copy(cepa_reads_core_s, cepa_reads_core_s_bk)
                
            in_fasta_r1 = cepa_reads_core_r1
            in_fasta_r2 = cepa_reads_core_r2
            shutil.copy(cepa_reads_core_r1, cepa_reads_core_r1_bk)
            shutil.copy(cepa_reads_core_r2, cepa_reads_core_r2_bk)
            
        for x in os.listdir(output):
            if x.endswith(".bk"):
                os.remove(x)
            if x.endswith(".bam"):
                os.remove(x)
        execute(f'gzip {cepa_reads_core_r1}')
        execute(f'gzip {cepa_reads_core_r2}')
        if in_fasta_s:
            execute(f'gzip {cepa_reads_core_s}')


if __name__ == '__main__':
    import argparse
    import os
    from glob import glob
    from SNDG import mkdir

    parser = argparse.ArgumentParser(description='Mapping to variant calls pipeline.')
    required = parser.add_argument_group('required arguments')
    required.add_argument('-r', '--ref_dirs', required=True)
    required.add_argument('-o', '--output', required=True)

    required.add_argument('-s', '--sample_name', required=True)
    required.add_argument('-rs', '--in_fasta_s', default=None)
    required.add_argument('-r2', '--in_fasta_r2', required=True)
    required.add_argument('-r1', '--in_fasta_r1', required=True)
    required.add_argument('-t', '--tmp_dir', default="/tmp")
    required.add_argument("--cpus", default=1)

    args = parser.parse_args()

    if not os.path.exists(args.ref_dirs):
        raise FileNotFoundError(f"{args.ref_dirs} does not exists")
    if args.in_fasta_s and not os.path.exists(args.in_fasta_s):
        raise FileNotFoundError(f"{args.in_fasta_s} does not exists")
    if not os.path.exists(args.in_fasta_r2):
        raise FileNotFoundError(f"{args.in_fasta_r2} does not exists")
    if not os.path.exists(args.in_fasta_r1):
        raise FileNotFoundError(f"{args.in_fasta_r1} does not exists")

    refs = glob(f'{args.ref_dirs}/*.fasta') + glob(f'{args.ref_dirs}/*.fna')

    if not refs:
        raise FileNotFoundError(f"no references detected in {args.ref_dirs}")

    mkdir(args.output)
    if not os.path.exists(args.output):
        raise FileNotFoundError(f"could not create {args.output}")
    
    if not args.output.endswith("/"):
        args.output = args.output + "/"
    
    mrca = MultiRefCoreAln(base_refs=refs, cpus=args.cpus)
    mrca.core_reads(args.output, args.sample_name,
                    args.in_fasta_r1, args.in_fasta_r2, args.in_fasta_s, args.tmp_dir)

