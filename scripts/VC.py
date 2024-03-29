#!/usr/bin/env python3

from glob import glob
import os
import subprocess as sp
import sys
from SNDG.Comparative.VariantSet import VariantSetUtils


if __name__ == '__main__':
    import argparse

    from tqdm import tqdm

    parser = argparse.ArgumentParser(description='Utils over genebank file')

    subparsers = parser.add_subparsers(help='commands', description='valid subcommands', dest='command', required=True)

    cmd = subparsers.add_parser('haplotypecall', help='creates raw haplotype call from a list of bams')
    cmd.add_argument('bam_folder')
    cmd.add_argument('vcfsraw_folder')
    cmd.add_argument('reference_file')
    cmd.add_argument("--gatk_env", default=os.environ.get("GATK_ENV",""))


    #cmd.add_argument('--new_lt', default=None, help="new locus tag. Use ONLY if no locus tag is present")

    cmd = subparsers.add_parser('combineraw', help='combine raw vcfs from a folder')
    cmd.add_argument('vcfsraw_folder')
    cmd.add_argument('combined_raw_vcf')
    cmd.add_argument('reference_file')
    cmd.add_argument("--gatk_env", default=os.environ.get("GATK_ENV",""))

    cmd = subparsers.add_parser('genotype_call', help='genotype call from raw vcf')
    cmd.add_argument('combined_raw_vcf')
    cmd.add_argument('combined_vcf')
    cmd.add_argument('reference_file')
    cmd.add_argument("--gatk_env", default=os.environ.get("GATK_ENV",""))

    cmd = subparsers.add_parser('render_tree', help='evaluates the best parameters for the phylogeny')
    cmd.add_argument('tree_newick')
    cmd.add_argument('tree_png')

    args = parser.parse_args()

    if args.command == "haplotypecall":
        if args.gatk_env:
            sys.stderr.write(f"GATK env: {args.gatk_env}")
        with tqdm(glob(f'{args.bam_folder}/*.bam')) as pbar:
            for bam in pbar:
                pbar.set_description(bam)
                vcf = os.path.basename(bam).replace(".bam",".gvcf.gz")
                cmd = f'''{args.gatk_env} gatk  HaplotypeCaller -ERC GVCF -R "{args.reference_file}" \
                          -ploidy 2 -I {bam} --output-mode EMIT_ALL_CONFIDENT_SITES \
                          -O "{args.vcfsraw_folder}/{vcf}" '''
                sp.call(cmd,shell=True)

    if args.command == "combineraw":
        raw_vcfs = [f'"{x}"' for x in glob(f'{args.vcfsraw_folder}/*.gvcf.gz')]
        vcfs = " --variant ".join(raw_vcfs)
        cmd = f'''{args.gatk_env} gatk CombineGVCFs -R "{args.reference_file}" -O "{args.combined_raw_vcf}" --variant {vcfs}'''
        sp.call(cmd,shell=True)

    if args.command == "genotype_call":
        cmd = f'''{args.gatk_env} gatk GenotypeGVCFs -R "{args.reference_file}" -ploidy 2 \
                    -V "{args.combined_raw_vcf}" -O "{args.combined_vcf}"'''
        sp.call(cmd,shell=True)
