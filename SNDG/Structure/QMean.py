import matplotlib
import shutil

matplotlib.use('Agg')

import tempfile

from ost.io import LoadPDB
from qmean import AssessModelQuality


class QMean:

    @staticmethod
    def assesment(pdb_path, output_dir=None):
        if not output_dir:
            output_dir = tempfile.mkdtemp(suffix="_qmean")
        pdb = LoadPDB(pdb_path)
        assessment = AssessModelQuality(pdb, output_dir=output_dir)
        shutil.rmtree(output_dir)
        return assessment[0].qmean4.__dict__


if __name__ == '__main__':
    from SNDG import init_log

    init_log()
    #assessment = QMean.assesment("test/P9WKI1_1gr0A.pdb")
    assessment = QMean.assesment("/data/databases/pdb/pdb/divided/ok/pdb4oke.ent")
    print(assessment)
