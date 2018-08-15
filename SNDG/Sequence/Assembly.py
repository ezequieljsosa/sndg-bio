"""

"""

from SNDG import execute, mkdir


class Assembly:

    @staticmethod
    def assemble_pe(libraries, out, trusted_contigs=None, untrusted_contigs=None, cov_cutoff=5):
        """

        :param libraries: list of tuples [(l1_r1_path,l1_r2_path,l1_single1_path,l1_single2_path,...),(l2_r1_path,...),... ]
        :param out: output dir
        :param trusted_contigs: trusted contigs path
        :param untrusted_contigs: untrusted contigs path
        :param cov_cutoff:
        :return:
        """
        assert libraries, "libraries cannot be empty"
        mkdir(out)

        template = """spades.py {libs} {tcont} {utcont} --cov-cutoff {cov_cutoff} -o {out}"""
        libs = ""
        for i, lib in enumerate(libraries, 1):
            libs += ' --pe{i}-1 "{r1}" --pe{i}-2 "{r2}" '.format(i=i, r1=lib[0], r2=lib[1])
            if len(lib) > 2:
                libs += " ".join([' --pe{i}-s "{ss}" '.format(i=i, ss=x) for x in lib[2:]])
        tcont = " --trusted-contigs " + trusted_contigs if trusted_contigs else ""
        utcont = " --untrusted-contigs " + untrusted_contigs if untrusted_contigs else ""
        execute(template, libs=libs, tcont=tcont, utcont=utcont, cov_cutoff=cov_cutoff, out=out)

    @staticmethod
    def quast(glob_exp, out, ref=None):
        execute("quast " + glob_exp + " -o " + out + ((" -R  " + ref) if ref else ""))

    def improve_assembly(self):
        pass


if __name__ == '__main__':
    import argparse
    import os
    from tqdm import tqdm
    from SNDG import init_log
    import Bio.SeqIO as bpio

    init_log()

    # parser = argparse.ArgumentParser(description='Assembly pipeline.')
    # required = parser.add_argument_group('required arguments')
    # required.add_argument('-o', '--work_dir', action='store', dest='work_dir', required=True)
    # required.add_argument('-S', '--strain', action='store', dest='strain', default="sample")
    # required.add_argument('-R1', '--read1', action='store', dest='read1', required=True)
    # required.add_argument('-R2', '--read2', action='store', dest='read2', required=True)
    #
    # args = parser.parse_args()
    # mkdir(args.work_dir)
    sm = "source /opt/PAGIT/sourceme.pagit"

    ref = "/data/projects/23staphylo/processed/refn315/genomic.fasta"
    refann = "/data/projects/23staphylo/processed/refn315/genomic.gb"
    for strain in tqdm(os.listdir("/data/projects/23staphylo/processed/mapping-n315")):
        base = "/data/projects/23staphylo/processed/denovo_1_improved/" + strain + "/"
        wd = base + "/runABACAS"
        contigs = "/data/projects/23staphylo/processed/denovo_1/" + strain + ".fasta"

        out_abacas = wd + "/" + os.path.basename(contigs) + "_" + os.path.basename(ref) + ".fasta"
        if not os.path.exists(out_abacas) or os.path.getsize(out_abacas) < 10:
            mkdir(wd)
            abacas = "perl $PAGIT_HOME/ABACAS/abacas.pl -r {ref} -q {contigs} -p nucmer -d > out.abacas.txt"
            execute("bash -c '{sm} ; " + abacas + "'", wd, ref=ref, contigs=contigs, sm=sm)

        assert os.path.exists(out_abacas), out_abacas

        wd = base + "/runIMAGE"
        mkdir(wd)
        image_out = base + "/contigs.fasta"
        if not os.path.exists(wd + "/ite6/new.fa") or os.path.getsize(wd + "/ite6/new.fa") < 10:
            r1 = "/data/projects/23staphylo/processed/trimmed2/" + strain + "_R1.fq.gz"
            r2 = "/data/projects/23staphylo/processed/trimmed2/" + strain + "_R2.fq.gz"

            execute("zcat {r1} > r_1.fastq", wd, r1=r1)
            execute("zcat {r2} > r_2.fastq", wd, r2=r2)
            image1 = "$PAGIT_HOME/IMAGE/image.pl -scaffolds {out_abacas} -prefix  r -iteration 1 -all_iteration 2 -dir_prefix ite -kmer 61 > out.image.txt"
            execute("bash -c '{sm};" + image1 + "'", wd, ref=ref, contigs=contigs, sm=sm, out_abacas=out_abacas)
            image2 = "$PAGIT_HOME/IMAGE/restartIMAGE.pl  ite2 49 2 partitioned > out.image2.txt"
            execute("bash -c '{sm};" + image2 + "'", wd, sm=sm)
            image3 = "$PAGIT_HOME/IMAGE/restartIMAGE.pl  ite4 41 2 partitioned > out.image3.txt"
            execute("bash -c '{sm};" + image3 + "'", wd, sm=sm)
            execute("rm *.fastq", wd)
            image4 = "contigs2scaffolds.pl ite6/new.fa ite6/new.read.placed 300 10 Res.image >> ./out.image.pl"
            execute("bash -c '{sm};" + image4 + "'", wd, ref=ref, contigs=contigs, sm=sm)
            with open(image_out,"w") as h:
                for r in bpio.parse(wd + "/ite6/new.fa", "fasta"):
                    r.id = strain + r.id.split(".")[0]
                    r.name = r.id
                    for f in r.features:
                        fe.qualifiers["note"] = fe.qualifiers["locus_tag"]

                    bpio.write(r, h , "fasta")
            execute("mv Res.image.* ../", wd)
            execute("rm -rf " + wd)

        if not os.path.exists(wd + "/out.ratt.txt") or os.path.getsize(wd + "/out.ratt.txt") < 10:
            wd = base + "/runRATT"
            mkdir(wd + "/embl")
            bpio.write(bpio.parse(refann, "gb"), wd + "/embl/ann.embl", "embl")
            ratt = "start.ratt.sh embl {image_out} Transfer Species > out.ratt.txt"
            execute("bash -c '{sm};" + ratt + "'", wd, image_out=image_out, sm=sm)
            #sed  -i 's/; ; ; ; ;/; SV 1; ; DNA; ; ;/' *.final.embl
            """"""

    """
    sed  -i 's/; ; ; ; ;/; SV 1; ; DNA; ; ;/' *.final.embl
    with open("./ann.gb","w") as h:
         for x in glob.glob("./runRATT/*.final.embl"):
             print x
             f = list(bpio.parse(x,"embl"))        
             for s in f:
                 s.id = s.id.split(".")[1]
                 s.name = s.id
            bpio.write(f,h,"gb")   
    
    """
