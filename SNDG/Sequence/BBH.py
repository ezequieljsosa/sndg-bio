"""

"""
import logging
from tqdm import tqdm
from collections import defaultdict
import Bio.SearchIO as bpsio
from SNDG.Sequence import read_blast_table

_log = logging.getLogger(__name__)


class BBH:

    @staticmethod
    def bbhs_from_blast(path_file1, path_file2, ident_threshold=80):
        """
        :param path_file1: blast result table format (6)
        :param path_file2: blast result table format (6)
        :return: list of bbh tuples
        """
        query_dict = defaultdict(lambda: {})
        for _, r in read_blast_table(path_file1).iterrows():

            if r.identity > ident_threshold:
                query_dict[r["query"]][r["hit"]] = 1

        bbhs = []
        for _, r in read_blast_table(path_file2).iterrows():
            if r.identity > ident_threshold:
                if r["query"] in query_dict[r["hit"]]:
                    bbhs.append((r["hit"], r["query"]))
        return bbhs

    @staticmethod
    def bbhs_from_lists(path_file1, path_file2, strict=False):
        """
        :param path_file1: blast result table format (6)
        :param path_file2: blast result table format (6)
        :return: list of bbh tuples
        """
        query_dict = defaultdict(lambda: {})
        for line in open(path_file1).readlines():
            try:
                x, y = line.strip().split("\t")[0:2]
                query_dict[x][y] = 1
            except:
                if strict:
                    raise Exception("parse Error: " + line)
                else:
                    _log.warn("parse Error: " + line)

        bbhs = []
        for line in open(path_file2).readlines():
            x, y = line.strip().split("\t")
            if x in query_dict[y]:
                bbhs.append([y, x])
        return bbhs


if __name__ == '__main__':
    import itertools
    import os

    # os.chdir("/mnt/Data/data/organismos/ILEX_PARA/blast")
    for x, y in itertools.product(["Genome", "Arabidopsis"], ["DEG", "KOG", "BUSCO"]):
        file_path1 = x + "_" + y + "2.blastout"
        file_path2 = y + "_" + x + "2.blastout"
        print (x, y)
        data = BBH.bbhs_from_blast(file_path1, file_path2, ident_threshold=0.8)
        data = set([x for x, y in data])
        print  len(data)
