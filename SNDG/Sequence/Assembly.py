"""

"""
import os
import shutil
from SNDG import execute, mkdir
import Bio.SeqIO as bpio
import multiprocessing

class Assembly:
    SPADES_DOCKER_IMAGE = 'ezequieljsosa/spades'

    @staticmethod
    def assemble_pe(r1: str, r2: str, out: str, name: str , ss: str = None, trusted_contigs: str = None,
                    untrusted_contigs: str = None, cov_cutoff: int = 5,tmp_dir:str="/tmp/",cpus=multiprocessing.cpu_count()):
        """

        :param out: output dir
        :param trusted_contigs: trusted contigs path
        :param untrusted_contigs: untrusted contigs path
        :param cov_cutoff:
        :return:
        """
        if tmp_dir == "/tmp/":
            tmp_dir = tmp_dir + name


        workdir1 = tmp_dir
        workdir2 = os.path.dirname(r1)
        assert workdir2 == os.path.dirname(r2), "r1 and r2 must be in the same directory"
        mkdir(out)

        mappings = f" -v {workdir1}:/out "
        in_dir = "/out/"
        if workdir1 != workdir2:
            mappings = mappings + f" -v {workdir2}:/in "
            in_dir = "/in/"

        template = """docker run  -u $(id -u):$(id -g) --rm -w /out {mappings} {image} spades.py \
                        {libs} {tcont} {utcont} -t {cpus} --isolate --cov-cutoff {cov_cutoff} -o /out """
        libs = ""

        i = 1
        r1_img = in_dir + r1.split(workdir2)[1]
        r2_img = in_dir + r2.split(workdir2)[1]
        libs += f' --pe{i}-1 "{r1_img}" --pe{i}-2 "{r2_img}" '
        if ss:
            ss_img = in_dir + ss.split(workdir2)[1]
            libs += f' --pe{i}-s "{ss_img}" '

        tcont = " --trusted-contigs " + trusted_contigs if trusted_contigs else ""
        utcont = " --untrusted-contigs " + untrusted_contigs if untrusted_contigs else ""
        cmd = template.format( libs=libs, tcont=tcont, utcont=utcont, cov_cutoff=cov_cutoff, out=out,
            mappings=mappings, image=Assembly.SPADES_DOCKER_IMAGE,cpus=cpus)
        print(cmd)
        # execute(cmd)
        # with open(f"{out}/scaffolds.fasta","w") as h:
        #     for r in bpio.parse(f"{tmp_dir}/scaffolds.fasta","fasta"):
        #         data = r.id.split("_") #"NODE_3_length_237403_cov_243.207"
        #         r.id = "|".join( [name + "_" + data[1],data[3],data[5].split(".")[0]] )
        #
        #         bpio.write(r,h,"fasta")
        # shutil.rmtree(f'"{tmp_dir}/{name}"')

    @staticmethod
    def quast(glob_exp, out, ref=None):
        execute("quast " + glob_exp + " -o " + out + ((" -R  " + ref) if ref else ""))

    def improve_assembly(self):
        pass


if __name__ == '__main__':
    import argparse
    from tqdm import tqdm
    assembly = Assembly()

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

