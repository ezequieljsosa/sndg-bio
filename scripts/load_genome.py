import argparse
import gzip
import Bio.SeqIO as bpio
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from tqdm import tqdm
import sys


parser = argparse.ArgumentParser()
parser.add_argument('-g', help='foo help',required=True)
args = parser.parse_args()


from BioSQL import BioSeqDatabase


server = BioSeqDatabase.open_database(driver="MySQLdb", user="root",
                                      passwd="docker", host="db", db="bioseqdb")

acc = "GCF_" +  args.g.split("_")[1]

if acc in server:
    sys.exit(0)

cs = []
with  gzip.open(args.g,"rt") as h:
    for x in tqdm(bpio.parse( h ,"gb")):
        cs.append(x)
assert (cs)

db = server.new_database(acc, description="")
server.commit()

count = db.load(tqdm(cs))
print (count)
server.commit()

db = server.new_database(acc  + "_prots", description="")
server.commit()
prots = []
iterator = iter(tqdm(cs))
while True:
    try:
        x = next(iterator)
    except StopIteration:
        break
    except Exception as e:
        print (e)
        continue
    for f in tqdm(x.features):
        if f.type == "CDS" and "translation" in f.qualifiers and "locus_tag" in f.qualifiers :
            lt = f.qualifiers["locus_tag"][0]
            record = SeqRecord(id=lt,name=lt,description="",seq=Seq(f.qualifiers["translation"][0]))
            prots.append(record)
        if len(prots) == 1000:
            count = db.load(prots)
            server.commit()
            prots = []
if prots:
    count = db.load(prots)
    server.commit()

"""
from glob import glob 
import subprocess as sp
from tqdm import tqdm 
import multiprocessing                
import time                           
iterator =   tqdm(glob("/out/*/*.gz"))
for x in iterator:             
    iterator.set_description(x.split("/")[-1])                                    
    def task():                                                    
        try:                                                       
            sp.call("python load_genome.py -g '" + x + "'", shell=True)
        except Exception as ex:
            print([x,ex])
    p = multiprocessing.Process(target=task)
    p.start()    
    p.join(5*60)
    if p.is_alive():
        p.terminate()
        p.join()    
"""