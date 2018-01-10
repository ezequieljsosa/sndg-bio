'''
Created on Oct 19, 2017

@author: eze
'''

import pandas as pd
import subprocess as sp
from Bio.SeqFeature import SeqFeature, FeatureLocation

def default_finfo(f):
    return f.qualifiers["NOTE"][0] if "NOTE" in f.qualifiers else (f.qualifiers["product"][0] if "product" in f.qualifiers else "")
def default_fname(f):
    return f.qualifiers["gene"][0] if "gene" in f.qualifiers else f.id

class CoverageAnalysis(object):
    '''
    samtools depth
    '''
    
    def __init__(self, depth_path=None):
        self.min_depth = 1 
        if depth_path:
            self.load_depth(depth_path)
    
    def run_sdepth(self,aln):
        depth_path = "/tmp/depth.txt"
        sp.call("samtools depth " + aln + " > " + depth_path,shell=True)
        self.load_depth(depth_path)
    
    def load_depth(self,depth_path):
        self.depth_df = pd.read_table(depth_path, names=["contig", "pos", "depth"])
    
    
    
    
    def feature_coverage(self, features, fname=default_fname, finfo=default_finfo,skip=["source","contig"]):
        
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
                { "type":f.type, "start":f.location.start,
                  "end":f.location.end, "strand":f.location.strand,
                "name":fname(f), "info":finfo(f), "hcov": round(len(data) * 1.0 / total,2) ,
                "depth_mean":toFloat(data.depth.mean()), "depth_min": toFloat(data.depth.min()),
                "depth_max":toFloat(data.depth.max())}, ignore_index=True)
    
    def feature_coverage_table(self,cfile="coverage.csv"):
        self.feature_cov_df.to_csv(cfile, columns=["type", "start", "end", "strand", "name", "info","hcov",
                                                    "depth_mean", "depth_min", "depth_max"],  index=False)
    
    def not_cov_features(self,cov_threshold,cov_fn=lambda r,ct:r.hcov >= ct):
        for _,r in self.feature_cov_df.iterrows():
            if not cov_fn(r,cov_threshold):
                yield SeqFeature(FeatureLocation(start=int(r["start"]),end=int(r["end"]),strand=int(r["strand"])), type=r["type"],  id=r["name"],
                                 qualifiers={"info":r["info"],"hcov":r["hcov"],"depth_mean":r["depth_mean"]
                                             ,"depth_min":r["depth_min"],"depth_max":r["depth_max"]} )
    
    def covered_features(self,cov_threshold,cfeatures="cfeatures.txt",ncfeatures="ncfeatures.txt",cov_fn=lambda r,ct:r.hcov >= ct ):
        with open(cfeatures,"w") as h:
            with open(ncfeatures,"w") as h2:
                for _,r in self.feature_cov_df.iterrows():
                    if cov_fn(r,cov_threshold):
                        h.write( r["name"]  + "\n")
                    else:
                        h2.write( r["name"] + "\n" )
            
    def gb_cov(self,annotation,fname=default_fname, finfo=default_finfo):
        import Bio.SeqIO as bpio
        features = next(bpio.parse(annotation,"gb")).features
        self.feature_coverage(features, fname=fname, finfo=finfo)
        
        