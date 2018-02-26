"""

"""
import gzip
import os
from glob import glob

import numpy as np
import pandas as pd
from tqdm import tqdm

import Bio.SeqIO as bpio
from SNDG import execute, mkdir


class FastQ:

    @staticmethod
    def rawstats(workdir, strains="", ref_size=None):
        columns = ["entry", "size_mb", "read_size_bp", "readcount"]
        if ref_size:
            columns.append("depth")
        dicts = []
        for strain in tqdm(strains):
            for filename in glob(workdir + "/" + strain + "*.gz"):
                entry = os.path.basename(filename).replace(".fastq.gz", "").replace(".fq.gz", "")
                total_bp = 0
                read_count = 0
                with gzip.open(filename) as h:
                    for read in bpio.parse(h, "fastq"):
                        read_count += 1
                        total_bp += len(read)

                size = os.path.getsize(filename) * 1.0 / 1024 / 1024
                if ref_size:
                    cov = total_bp / ref_size
                    dicts.append((entry, size, total_bp / read_count, read_count, cov,))
                else:
                    dicts.append((entry, size, total_bp / read_count, read_count,))

        df = pd.DataFrame.from_records(dicts, index=None, columns=columns)
        return df

    @staticmethod
    def filteredstats(workdir, strains="", ref_size=None, rawdf=None):
        columns = ["entry", "f_size_mb", "f_total_bp", "f_readcount"]

        if ref_size:
            columns.append("f_depth")
        dicts = []

        for strain in tqdm(strains):
            for filename in glob(workdir + "/" + strain + "*.gz"):

                entry = os.path.basename(filename).replace(".fastq.gz", "").replace(".fq.gz", "")
                total_bp = 0
                read_count = 0
                with gzip.open(filename) as h:
                    for read in bpio.parse(h, "fastq"):
                        read_count += 1
                        total_bp += len(read)

                size = os.path.getsize(filename) * 1.0 / 1024 / 1024

                if ref_size:
                    cov = total_bp * 1.0 / ref_size
                    dicts.append((entry, size, total_bp, read_count, cov,))
                else:
                    dicts.append((entry, size, total_bp, read_count,))

        df = pd.DataFrame.from_records(dicts, index=None, columns=columns)

        try:
            return pd.merge(rawdf, df, on=["entry"], how='right')
        except:
            return df

    @staticmethod
    def sumarize_stats_df(strains, df):
        """
        Takes a df with the reads file name in "entry" column, and colapses it
        to on row per sample. It assumes that sample name is at the beginning of the read file name.
        :param df:
        :return:
        """
        srows = []

        df2 = df.copy()
        df2["strain"] = [[y for y in strains if y in x][0] for x in df2.entry]

        for strain, group in df2.groupby(["strain"]):

            rows = [x for _, x in group.iterrows()]

            new_row = [strain]
            for c, fn in [("read_size_bp", lambda xs: sum(xs) / len(xs)),
                          ("f_total_bp", sum), ("depth", sum),
                          ("f_depth", sum), ("readcount", sum), ("f_readcount", sum)]:
                value =  round(fn([x[c] for x in rows if not np.isnan(x[c])]),2)
                new_row.append(value)
            srows.append(new_row)
        cols = ["strain", "read_size_bp", "f_total_bp",
                "depth", "f_depth", "readcount", "f_readcount"]
        return pd.DataFrame.from_records(srows, columns=cols)

    @staticmethod
    def fastqc(source_dir, dst_dir):
        mkdir(dst_dir)
        for filename in tqdm(sorted(glob(source_dir + "/*"))):
            execute("fastqc  {src} -q --extract -o {dst}",
                    src=filename, dst=dst_dir)

    @staticmethod
    def trim(strains, source_dir, dst_dir, clip="", headcrop=13, quality=20, windowsize=4, minlen=36):
        """

        :param strains:
        :param source_dir:
        :param dst_dir:
        :param clip: ILLUMINACLIP:../data/external/NexteraPE-PE.fa:2:30:10
        :param headcrop:
        :param quality:
        :param windowsize:
        :param minlen:
        :return:
        """
        mkdir(dst_dir)
        with tqdm(strains) as pbar:
            for strain in pbar:
                filenames = [os.path.basename(x) for x in glob(source_dir + "/" + strain + "*.gz")]
                cmd = """trimmomatic PE \
    ../data/raw/{r1} ../data/raw/{r2}\
    {dst}/{r1} {dst}/{r1u} \
    {dst}/{r2} {dst}/{r2u} \
     {clip} HEADCROP:{headcrop} \
     LEADING:{quality} TRAILING:{quality} SLIDINGWINDOW:{windowsize}:20 MINLEN:{minlen}"""
                cmd = cmd.format(clip=clip, quality=quality, windowsize=windowsize, minlen=minlen,
                                 headcrop=headcrop, dst=dst_dir, r1=filenames[0], r2=filenames[1],
                                 r1u=filenames[0].replace(".fastq.gz", "_unpaired.fastq.gz"),
                                 r2u=filenames[1].replace(".fastq.gz", "_unpaired.fastq.gz"))
                pbar.set_description(cmd)
                execute(cmd)

                execute(" gunzip " + dst_dir + "/" + filenames[0].replace(".fastq.gz", "_unpaired.fastq.gz"))
                execute(" gunzip " + dst_dir + "/" + filenames[1].replace(".fastq.gz", "_unpaired.fastq.gz"))
                ur1 = dst_dir + "/" + filenames[0].replace(".fastq.gz", "_unpaired.fastq")
                ur2 = dst_dir + "/" + filenames[1].replace(".fastq.gz", "_unpaired.fastq")
                execute("cat " + ur1 +
                        " " + ur2 +
                        " > " + dst_dir + "/" + strain + "_unpaired.fastq")
                os.remove(ur1)
                os.remove(ur2)
                execute("gzip " + dst_dir + "/" + strain + "_unpaired.fastq")
