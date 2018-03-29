from argparse import ArgumentParser, RawDescriptionHelpFormatter

import matplotlib
import shutil

import os
from Bio.PDB import PDBParser, is_aa
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import seq1

matplotlib.use('Agg')

from SNDG import execute, mkdir

import tempfile
import Bio.SeqIO as bpio
from ost.io import LoadPDB
from qmean import AssessModelQuality
from qmean import PSIPREDHandler
from qmean import ACCPROHandler
import multiprocessing


class QMean:

    @staticmethod
    def psipred(fasta, path="./", cmd="/opt/psipred/runpsipred",cpus=multiprocessing.cpu_count()):
        """
        #runpsipred example.fasta --> example.horiz
        # PSIPRED HFORMAT (PSIPRED V3.5)

Conf: 928999937998289999999999961696258972566893341566778999832667
Pred: CEEEEEEECCCCCHHHHHHHHHHHHHHCCCEEEEEECCCCCCCCCCCCCCEEEEEECCCC
  AA: PKALIVYGSTTGNTEYTAETIARQLANAGYEVDSRDAASVEAGGLFEGFDLVLLGCSTWG
              10        20        30        40        50        60

        :param fasta:
        :return:
        """
        execute(cmd + " " + fasta + " " + str(cpus), wd=path)
        horiz = fasta.replace(".fasta", ".horiz")
        pred = ""
        conf = ""
        with open(horiz) as h:
            for x in h.readlines():
                if x.startswith("Pred:"):
                    pred += x.split(" ")[1].strip()
                if x.startswith("Conf:"):
                    conf += x.split(" ")[1].strip()

        return (pred, conf)

    @staticmethod
    def accpro(fasta, path, cmd="/opt/sspro4/bin/predict_acc.sh"):
        """
        #../bin/predict_ssa.sh 1aqta.fasta 1aqta.test
        cat /opt/sspro4/test/1aqta.acc6

        1aqta_fastaalg
STYHLDVVSAEQQMFSGLVEKIQVTGSEGELGIYPGHAPLLTAIKPGMIRIVKQHGHEEFIYLSGGILEVQPGNVTVLADTAIRGQDLDEARAMEAKRKAEEHISSSHGDVDYAQASAELAKAIAQLRVIELTKK
eebebbbbbbeeebbeeebeebbbebeebbbbbbbebbbbbbbbebbbbbbebeeeeebbbbbbbbbbbbeeeebbbbbbbbeeeeebeeeebeebbeebbeebeeeeeeeebeebeebbeebbebbebbeeeee


        :param fasta:
        :param path:
        :param cmd:
        :return:
        """
        acc6 = fasta.replace(".fasta", ".acc6")
        if not os.path.exists(acc6):
            execute(cmd + " " + fasta + " " + acc6, wd=path)

        with open(acc6) as h:
            return h.readlines()[2].strip()

    @staticmethod
    def assesment(pdb_path, fasta_path=None, output_dir=None,
                  accpro_path=None, psipred_path=None,cpus=multiprocessing.cpu_count()):

        if not output_dir:
            output_dir = tempfile.mkdtemp(suffix="_qmean")

        data = dict()
        if not fasta_path and (accpro_path or psipred_path):
            fasta_path = output_dir + "/seq.fasta"

            p = PDBParser(PERMISSIVE=True, QUIET=True)
            seq = "".join([seq1(residue.resname) for residue in p.get_structure("x", pdb_path).get_residues()
                           ])  # if is_aa(residue)
            with open(fasta_path, "w") as h:
                h.write(">seq\n")
                h.write(seq)

            data["seq"] = seq

        psipred_handler = None
        accpro_handler = None

        if accpro_path:
            data["acc"] = QMean.accpro(fasta_path, output_dir, accpro_path)
            accpro_handler = ACCPROHandler(data)

        if psipred_path:
            data["ss"], data["conf"] = QMean.psipred(fasta_path, output_dir, psipred_path,cpus)
            psipred_handler = PSIPREDHandler(data)

        pdb = LoadPDB(pdb_path)

        if psipred_handler and accpro_handler:
            assessment = AssessModelQuality(pdb, output_dir=output_dir,
                                            psipred=psipred_handler, accpro=accpro_handler)
        elif psipred_handler and not accpro_handler:
            assessment = AssessModelQuality(pdb, output_dir=output_dir,
                                            psipred=psipred_handler)
        elif not psipred_handler and accpro_handler:
            assessment = AssessModelQuality(pdb, output_dir=output_dir,
                                            accpro=accpro_handler)
        else:
            assessment = AssessModelQuality(pdb, output_dir=output_dir)

        #shutil.rmtree(output_dir)
        result = {}
        for x in assessment[0].all_scores:
            result[x.name + "_norm"] = x.norm
            result[x.name + "_zscore"] = x.z_score
        result["residues"] = {}
        for row in assessment[1].score_table.rows:
            r = {f: row[i] for i, f in enumerate(assessment[1].score_table.col_names[4:], 4)}
            result["residues"][row[0] + "_" + str(row[2]) + "_" + str(row[3])] = r
        return result


if __name__ == '__main__':
    from SNDG import init_log

    init_log()

    parser = ArgumentParser(formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-acc", "--accpro", default="/opt/sspro4/bin/predict_acc.sh")
    parser.add_argument("-psi", "--psipred", default="/opt/psipred/runpsipred")
    parser.add_argument("-i", "--inputpdb", required=True)
    parser.add_argument("-o", "--outdir", default="./")
    parser.add_argument( "--cpus", default=multiprocessing.cpu_count())

    args = parser.parse_args()

    # "/data/databases/pdb/divided/ok/pdb4oke.ent"
    mkdir(args.outdir)
    assessment = QMean.assesment(args.inputpdb, output_dir=args.outdir,
                                 accpro_path=args.accpro, psipred_path=args.psipred,cpus=args.cpus)

    print(assessment)
