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

        template = """spades {libs} {tcont} {utcont} --cov-cutoff {cov_cutoff} -o {out}"""
        libs = ""
        for i, lib in enumerate(libraries, 1):
            libs += " --pe{i}-1 {r1} --pe{i}-2 {r2} ".format(i=i, r1=lib[0], r2=lib[1])
            if len(lib) > 2:
                libs += " ".join([" --pe{i}-s {ss} ".format(i=i, ss=x) for x in lib[2:]])
        tcont = " --trusted-contigs " + trusted_contigs if trusted_contigs else ""
        utcont = " --untrusted-contigs " + untrusted_contigs if untrusted_contigs else ""
        execute(template, libs=libs, tcont=tcont, utcont=utcont, cov_cutoff=cov_cutoff, out=out)

    @staticmethod
    def quast(glob_exp,out,ref = None):
        execute("quast " + glob_exp + " -o " + out + ((" -R  " + ref) if ref else "")  )
