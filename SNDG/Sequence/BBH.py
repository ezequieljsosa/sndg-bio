"""

"""
import logging
from tqdm import tqdm
from collections import defaultdict
import Bio.SearchIO as bpsio

_log = logging.getLogger(__name__)

class BBH:

    @staticmethod
    def bbhs_from_blast(path_file1, path_file2,ident_threshold=80):
        """
        :param path_file1: blast result table format (6)
        :param path_file2: blast result table format (6)
        :return: list of bbh tuples
        """
        query_dict = defaultdict(lambda :  {})
        for query in bpsio.parse(path_file1, "blast-tab"):
            for hit in query:
                hsp = hit[0]
                if hsp.ident_pct > ident_threshold:
                    query_dict[query.id][hit.id] = 1

        bbhs = []
        for query in bpsio.parse(path_file2, "blast-tab"):
            for hit in query:
                hsp = hit[0]
                if hsp.ident_pct > ident_threshold:
                    if query.id in query_dict[hit.id]:
                        bbhs.append(hit.id,query.id)
        return bbhs

    @staticmethod
    def bbhs_from_lists(path_file1, path_file2,strict=False):
        """
        :param path_file1: blast result table format (6)
        :param path_file2: blast result table format (6)
        :return: list of bbh tuples
        """
        query_dict = defaultdict(lambda :  {})
        for line in open(path_file1).readlines():
            try:
                x,y = line.strip().split("\t")[0:2]
                query_dict[x][y] = 1
            except:
                if strict:
                    raise Exception("parse Error: " + line)
                else:
                    _log.warn("parse Error: " + line)

        bbhs = []
        for line in open(path_file2).readlines():
            x,y = line.strip().split("\t")
            if x in query_dict[y]:
                bbhs.append([y,x])
        return bbhs

if __name__ == '__main__':
    import itertools
    import os
    os.chdir("/mnt/Data/data/organismos/ILEX_PARA/blast")
    for x,y in itertools.product( ["Arabidopsis","Transcriptome","Genome"],["DEG","KOG","BUSCO"] ) :
         file_path1 = x + "_" + y + ".blastout.hits"
         file_path2 = y + "_" + x + ".blastout.hits"
         print (x,y)
         print  len( BBH.bbhs_from_lists( file_path1,file_path2))
