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
        result = {}
        for x in assessment[0].all_scores:
            result[x.name + "_norm"] = x.norm
            result[x.name + "_zscore"] = x.z_score
        result["residues"] = {}
        for row in assessment[1].score_table.rows:
            r = {f:row[i] for i,f in enumerate(assessment[1].score_table.col_names[4:],4) }
            result["residues"][  row[0]  + "_" +  str(row[2]) + "_" +  str(row[3])] = r
        return result


if __name__ == '__main__':
    from SNDG import init_log

    init_log()
    #assessment = QMean.assesment("test/P9WKI1_1gr0A.pdb")
    assessment = QMean.assesment("/data/databases/pdb/pdb/divided/ok/pdb4oke.ent")
    print(assessment)
