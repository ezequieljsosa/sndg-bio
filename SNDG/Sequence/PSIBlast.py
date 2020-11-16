import Bio.SeqIO as bpio
from glob import glob
import argparse
import subprocess as sp


class PsiBlast:
    pass

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Create PSSMs for a fasta file')
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--db', required=True)
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--cpu', required=True)

    args = parser.parse_args()
    iterations = 3

    for seq in bpio.parse(args.fasta,"fasta"):
        output = args.outdir + "/" + seq.id + ".pssm"
        in_seq = f"/tmp/{seq.id}.fasta"
        bpio.write(seq,in_seq,"fasta")
        cmd = f"psiblast -query {in_seq} -db {args.db} -num_threads {args.cpu} -out_pssm {output} -evalue 0.0001 -num_iterations {iterations} -outfmt 6 > /dev/null"
        sp.run(cmd,shell=True)

