'''
Created on Oct 19, 2017

@author: eze
'''

import os
import tempfile

import pandas as pd

import Bio.SeqIO as bpio
from Bio.SeqFeature import SeqFeature, FeatureLocation

from SNDG import execute


def default_finfo(f):
    return f.qualifiers["NOTE"][0] if "NOTE" in f.qualifiers else (
        f.qualifiers["product"][0] if "product" in f.qualifiers else "")


def default_fname(f):
    return f.qualifiers["gene"][0] if "gene" in f.qualifiers else f.id


class CoverageAnalysis(object):
    '''
    samtools depth
    '''

    def __init__(self, depth_path="/tmp/depth.txt", min_depth=1):
        self.min_depth = min_depth
        self.depth_df = pd.DataFrame()
        self.depth_path = depth_path
        if os.path.exists(depth_path) and os.path.getsize(depth_path) > 100:
            self.load_depth()

    def run_sdepth(self, aln):
        execute("samtools depth -aa " + aln + " > " + self.depth_path)
        self.load_depth()

    def load_depth(self):
        self.depth_df = pd.read_table( self.depth_path, names=["contig", "pos", "depth"])

    def aligned_reads(self,aln):
        tmp = tempfile.mktemp()
        execute("samtools flagstat " + aln + " > "  + tmp)
        with open(tmp) as h:
            value = h.readline().split(" ")[0]
        return int(value)

    def horizontal_coverage(self):
        return  (1.0 * self.depth_df[self.depth_df.depth >= self.min_depth].size) / self.depth_df.size

    def feature_coverage(self, features, fname=default_fname, finfo=default_finfo, skip=["source", "contig"]):

        def toFloat(x):
            try:
                return float(x)
            except:
                return 0

        self.feature_cov_df = pd.DataFrame()
        for f in features:
            if f.type in skip:
                continue
            total = f.location.end - f.location.start
            data = self.depth_df[(self.depth_df.pos >= f.location.start)
                                 & (self.depth_df.pos < f.location.end)
                                 & (self.depth_df.depth > self.min_depth)]

            self.feature_cov_df = self.feature_cov_df.append(
                {"type": f.type, "start": f.location.start,
                 "end": f.location.end, "strand": f.location.strand,
                 "name": fname(f), "info": finfo(f), "hcov": round(len(data) * 1.0 / total, 2),
                 "depth_mean": toFloat(data.depth.mean()), "depth_min": toFloat(data.depth.min()),
                 "depth_max": toFloat(data.depth.max())}, ignore_index=True)

    def feature_coverage_table(self, cfile="coverage.csv"):
        self.feature_cov_df.to_csv(cfile, columns=["type", "start", "end", "strand", "name", "info", "hcov",
                                                   "depth_mean", "depth_min", "depth_max"], index=False)

    def not_cov_features(self, cov_threshold, cov_fn=lambda r, ct: r.hcov >= ct):
        for _, r in self.feature_cov_df.iterrows():
            if not cov_fn(r, cov_threshold):
                yield SeqFeature(FeatureLocation(start=int(r["start"]), end=int(r["end"]), strand=int(r["strand"])),
                                 type=r["type"], id=r["name"],
                                 qualifiers={"info": r["info"], "hcov": r["hcov"], "depth_mean": r["depth_mean"]
                                     , "depth_min": r["depth_min"], "depth_max": r["depth_max"]})

    def covered_features(self, cov_threshold, cfeatures="cfeatures.txt", ncfeatures="ncfeatures.txt",
                         cov_fn=lambda r, ct: r.hcov >= ct):
        with open(cfeatures, "w") as h:
            with open(ncfeatures, "w") as h2:
                for _, r in self.feature_cov_df.iterrows():
                    if cov_fn(r, cov_threshold):
                        h.write(r["name"] + "\n")
                    else:
                        h2.write(r["name"] + "\n")

    def gb_cov(self, annotation, fname=default_fname, finfo=default_finfo):
        features = next(bpio.parse(annotation, "gb")).features
        self.feature_coverage(features, fname=fname, finfo=finfo)
