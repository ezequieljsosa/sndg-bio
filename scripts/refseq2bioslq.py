import argparse
import gzip

import Bio.SeqIO as bpio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument('-refseq', help='refseq gzipped genbank file', required=True)
parser.add_argument('-db', help='database', default="bioseqdb")
parser.add_argument('-user', help='database user', default="root")
parser.add_argument('-host', help='database host', default="localhost")
parser.add_argument('-passwd', help='database password', required=True)
args = parser.parse_args()

from BioSQL import BioSeqDatabase

cs = []
with  gzip.open(args.refseq, "rt") as h:
    for x in tqdm(bpio.parse(h, "gb")):
        cs.append(x)
assert (cs)

server = BioSeqDatabase.open_database(driver="MySQLdb", user=args.user,
                                      passwd=args.passwd, host=args.host, db=args.db)

acc = "GCF_" + args.refseq.split("_")[1]
db = server.new_database(acc, description="")
server.commit()

count = db.load(tqdm(cs))
print (count)
server.commit()

db = server.new_database(acc + "_prots", description="")
server.commit()
prots = []
for x in tqdm(cs):
    for f in tqdm(x.features):
        if f.type == "CDS" and "translation" in f.qualifiers and "locus_tag" in f.qualifiers:
            lt = f.qualifiers["locus_tag"][0]
            record = SeqRecord(id=lt, name=lt, description="", seq=Seq(f.qualifiers["translation"][0]))
            prots.append(record)
count = db.load(tqdm(prots))
server.commit()
