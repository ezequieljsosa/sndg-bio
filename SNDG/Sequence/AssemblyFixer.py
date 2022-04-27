"""

"""
import os
import shutil
from SNDG import execute, mkdir
import Bio.SeqIO as bpio
import multiprocessing
from itertools import groupby
from SNDG.Sequence import smart_parse
from collections import  defaultdict

class BlastHit:

    def __init__(self, qseqid,qstart,qend,sseqid,sstart,send,pident,length,qcovhsp,slen):
        self.qseqid = qseqid
        self.qstart = qstart
        self.qend = qend
        self.sseqid = sseqid
        self.sstart = sstart
        self.send = send
        self.pident = pident
        self.length = length
        self.qcovhsp = qcovhsp
        self.slen = slen
        self.genes = []

    def addgene(self,gene):
        self.genes.append(gene)










class AssemblyFixer:

    @staticmethod
    def fromFile(file_path):
         return AssemblyFixer(bpio.to_dict(smart_parse(file_path)))

    def __init__(self, assembly):
        self.assembly = assembly


    def whole_protein_hint(self,trusted_protein):
        """
        The protein is splited between different contigs, but we trust it is on the assembly
        :return:
        """

    def primers_hint(self,trusted_protein):
        """
        we know that the region amplifies, and it is
        :return:
        """

    def scaffold_hint(self,hits,gene_order):
        """
        join the contigs using scaffold/reference hits
        """
        # sort hits


    def gene_order_hits(self,hits,gene_order,context=100):
        # sort hits
        hits_order = []
        contig_genes = defaultdict(list)
        for gene_data in gene_order:
            for hit in hits:
                if gene_data.id  == hit.qseqid:
                    hits_order.append(hit.sseqid)
                    hit.addgene(gene_data)

        hits_order_no_repeated = [x[0] for x in groupby(hits_order,key=lambda y:y.contig)]
        first = hits_order_no_repeated[0]
        gene_data = first.first_gene()
        seqs = [first.get_begining(gene_data)]

        for hit in hits_order_no_repeated[1:-1]:
            if hit.coverage > 80:
                gene_data = hit.first_gene()
                hit.sorted_contig(gene_data)

        gene_data = hits_order_no_repeated[-1].last_gene()
        seqs.append( first.get_end(gene_data) )

if __name__ == '__main__':
    import argparse
    from tqdm import tqdm
    assembly = AssemblyFixer()

    parser = argparse.ArgumentParser(description='Assembly pipeline.')
    required = parser.add_argument_group('required arguments')
    required.add_argument('-o', '--output',default=None)
    required.add_argument('-s', '--sample',  default="sample")
    required.add_argument('-r1', '--read1', required=True)
    required.add_argument('-r2', '--read2',  required=True)
    required.add_argument('-u', '--unpaired',  default=None)


    args = parser.parse_args()
    if not args.output:
        args.output = "./"  + args.sample
    assembly.assemble_pe(args.read1,args.read2,args.output,args.sample,ss=args.unpaired)

    # 3156, 3166, 3157,  3172

    # mkdir(args.work_dir)
    # sm = "source /opt/PAGIT/sourceme.pagit"
    #
    # ref = "/data/projects/23staphylo/processed/refn315/genomic.fasta"
    # refann = "/data/projects/23staphylo/processed/refn315/genomic.gb"
    # for strain in tqdm(os.listdir("/data/projects/23staphylo/processed/mapping-n315")):
    #     base = "/data/projects/23staphylo/processed/denovo_1_improved/" + strain + "/"
    #     wd = base + "/runABACAS"
    #     contigs = "/data/projects/23staphylo/processed/denovo_1/" + strain + ".fasta"
    #
    #     out_abacas = wd + "/" + os.path.basename(contigs) + "_" + os.path.basename(ref) + ".fasta"
    #     if not os.path.exists(out_abacas) or os.path.getsize(out_abacas) < 10:
    #         mkdir(wd)
    #         abacas = "perl $PAGIT_HOME/ABACAS/abacas.pl -r {ref} -q {contigs} -p nucmer -d > out.abacas.txt"
    #         execute("bash -c '{sm} ; " + abacas + "'", wd, ref=ref, contigs=contigs, sm=sm)
    #
    #     assert os.path.exists(out_abacas), out_abacas
    #
    #     wd = base + "/runIMAGE"
    #     mkdir(wd)
    #     image_out = base + "/contigs.fasta"
    #     if not os.path.exists(wd + "/ite6/new.fa") or os.path.getsize(wd + "/ite6/new.fa") < 10:
    #         r1 = "/data/projects/23staphylo/processed/trimmed2/" + strain + "_R1.fq.gz"
    #         r2 = "/data/projects/23staphylo/processed/trimmed2/" + strain + "_R2.fq.gz"
    #
    #         execute("zcat {r1} > r_1.fastq", wd, r1=r1)
    #         execute("zcat {r2} > r_2.fastq", wd, r2=r2)
    #         image1 = "$PAGIT_HOME/IMAGE/image.pl -scaffolds {out_abacas} -prefix  r -iteration 1 -all_iteration 2 -dir_prefix ite -kmer 61 > out.image.txt"
    #         execute("bash -c '{sm};" + image1 + "'", wd, ref=ref, contigs=contigs, sm=sm, out_abacas=out_abacas)
    #         image2 = "$PAGIT_HOME/IMAGE/restartIMAGE.pl  ite2 49 2 partitioned > out.image2.txt"
    #         execute("bash -c '{sm};" + image2 + "'", wd, sm=sm)
    #         image3 = "$PAGIT_HOME/IMAGE/restartIMAGE.pl  ite4 41 2 partitioned > out.image3.txt"
    #         execute("bash -c '{sm};" + image3 + "'", wd, sm=sm)
    #         execute("rm *.fastq", wd)
    #         image4 = "contigs2scaffolds.pl ite6/new.fa ite6/new.read.placed 300 10 Res.image >> ./out.image.pl"
    #         execute("bash -c '{sm};" + image4 + "'", wd, ref=ref, contigs=contigs, sm=sm)
    #         with open(image_out, "w") as h:
    #             for r in bpio.parse(wd + "/ite6/new.fa", "fasta"):
    #                 r.id = strain + r.id.split(".")[0]
    #                 r.name = r.id
    #                 for f in r.features:
    #                     fe.qualifiers["note"] = fe.qualifiers["locus_tag"]
    #
    #                 bpio.write(r, h, "fasta")
    #         execute("mv Res.image.* ../", wd)
    #         execute("rm -rf " + wd)
    #
    #     if not os.path.exists(wd + "/out.ratt.txt") or os.path.getsize(wd + "/out.ratt.txt") < 10:
    #         wd = base + "/runRATT"
    #         mkdir(wd + "/embl")
    #         bpio.write(bpio.parse(refann, "gb"), wd + "/embl/ann.embl", "embl")
    #         ratt = "start.ratt.sh embl {image_out} Transfer Species > out.ratt.txt"
    #         execute("bash -c '{sm};" + ratt + "'", wd, image_out=image_out, sm=sm)
    #         # sed  -i 's/; ; ; ; ;/; SV 1; ; DNA; ; ;/' *.final.embl
    #         """"""

