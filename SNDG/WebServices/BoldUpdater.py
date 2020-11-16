'''
Created on Feb 21, 2017

@author: eze
'''

import os
import json
import pymongo

from SNDG.WebServices import download_file

class Bold(object):
    '''
    classdocs
    '''

    def __init__(self, db,datadir):
        self.dbname = db
        self.datadir = datadir

    def download(self,datadir):
        self.bolddir = datadir + "bold/"
    
        if not os.path.exists(self.bolddir):
            os.makedirs(self.bolddir)     
        os.chdir(self.bolddir)    
        download_file(  "http://www.boldsystems.org/index.php/API_Public/combined?geo=Argentina&format=json" , "barcodes.json")
    
    def process(self): 
        self.bolddir = self.datadir + "bold/"
        db = pymongo.MongoClient()[self.dbname]
    
        if not os.path.exists(self.bolddir):
            os.makedirs(self.bolddir)     
        os.chdir(self.bolddir)    
        
        print (db.barcodes.remove({},multi=True))
         
        data = json.load(open("barcodes.json"))["bold_records"]["records"].values()
        for d in data:
            if "sequences" in d:
                for x in  ["species","genus","subfamily" ,"family","order"  ,"class" ,"phylum"]:
                        if x in d["taxonomy"]:
                            desc = d["taxonomy"][x]["taxon"]["name"]
                            tax = d["taxonomy"][x]["taxon"]["taxID"]
                            break   
                d["organism"] = desc
                
                d["description"] = (d["sequences"]["sequence"][0]["markercode"] + " " 
                                                           + d["specimen_identifiers"]["institution_storing"])
                
                d["tax"] = int(tax)
                db.barcodes.insert(d)
        
        

if __name__ == "__main__":
    Bold("saureus","/data/databases/").process() 
    print ("Ok")
        
        
        
       