'''
Created on Nov 8, 2017

@author: eze
'''
import os
import shutil

import Bio.SeqIO as bpio
import pymongo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import seq1
from SNDG import execute


class SNDGUpdater(object):
    '''
    
    '''

    def __init__(self, db):
        self.dbname = db

    def download(self, config, datadir):
        pass

    def process(self):

        db = pymongo.MongoClient()[self.dbname]
        pdbdb = pymongo.MongoClient()["pdb"]

        db.barcodes.remove({"sequences": {"$exists": False}})

        barcodes_path = "/data/databases/sndg/barcodes.fasta"
        with open(barcodes_path, "w") as h:
            for bc in db.barcodes.find():
                for x in ["species", "genus", "subfamily", "family", "order", "class", "phylum"]:
                    if x in bc["taxonomy"]:
                        desc = bc["taxonomy"][x]["taxon"]["name"]
                        break
                record = SeqRecord(description="", id=bc["processid"] + "||" + desc,
                                   seq=Seq(bc["sequences"]["sequence"][0]["nucleotides"].replace("-", "")))

                bpio.write(record, h, "fasta")
        execute("makeblastdb -dbtype nucl -in " + barcodes_path)

        seq_col_map = {str(x["_id"]): x["description"] if x["description"] else x["organism"]
                       for x in db.sequence_collection.find({}, {"description": 1, "organism": 1})}
        proteins_path = "/data/databases/sndg/proteins.fasta"
        with open(proteins_path, "w") as h:

            for p in db.proteins.find({}, {"seq_collection_id": 1, "seq": 1, "name": 1}):
                col_id = str(p["seq_collection_id"])
                if (col_id in seq_col_map) and (len(p["seq"]) > 10):
                    bpio.write(SeqRecord(
                        id=str(p["_id"]) + "||" + col_id + "||" + seq_col_map[col_id].replace(" ", "_") + "||" + p[
                            "name"],
                        description="", seq=Seq(p["seq"])), h, "fasta")
        execute("makeblastdb -dbtype prot -in " + proteins_path)

        structures_path = "/data/databases/sndg/structures.fasta"
        with open(structures_path, "w") as h:
            for p in pdbdb.structures.find({},
                                           {"chains.residues.compound": 1, "chains.name": 1, "organism": 1, "name": 1}):
                for c in p["chains"]:
                    col_id = p["name"] + "_" + c["name"]
                    seq = "".join([seq1(x["compound"]) for x in c["residues"]])
                    desc = ""
                    if "organism" in p:
                        desc = p["organism"].lower().replace(" ", "_")
                    if len(seq) > 10:
                        bpio.write(SeqRecord(id=str(p["_id"]) + "||" + col_id + "||" + desc,
                                             seq=Seq(seq)), h, "fasta")
        execute("makeblastdb -dbtype prot -in " + structures_path)

        # cnx = connector.connect(user='root', database='biosql',password='mito')
        # try:
        #     nucleotides_path = "/data/databases/sndg/nucleotides.fasta"
        #     with open(nucleotides_path,"w") as h:
        #             for contig in db.contig_collection.find({},{"name":1,"organism":1}):
        #                 cursor = cnx.cursor()
        #                 query = ("""SELECT DISTINCT be.identifier,bs.seq
        #                     FROM  bioentry be, biosequence bs
        #                     WHERE be.identifier LIKE '""" + contig["name"] + """%'
        #                         AND be.bioentry_id = bs.bioentry_id;""" )
        #
        #                 cursor.execute(query, ())
        #                 seqs = list(cursor)
        #                 for s in seqs:
        #                     bpio.write(SeqRecord(id=str(s[0].decode("utf-8"))  + "||" + contig["organism"] + "||" ,
        #                                      seq=Seq(s[1])), h, "fasta")
        #                 cursor.close()
        #     execute("makeblastdb -dbtype nucl -in " + nucleotides_path)
        # finally:
        #     cnx.close()

        print("ok")


if __name__ == "__main__":
    SNDGUpdater("saureus").process()

    stats_dir = "/data/xomeq/krona/"
    db = pymongo.MongoClient().saureus
    dbpdb = pymongo.MongoClient().pdb
    # cols = {"struct":"structures","prot":"proteins", "barcode":"barcodes"}

    os.chdir("/data/databases/sndg")

    with open("tax_tool.tsv", "w") as h:
        data = list(db.tools.aggregate([{"$group": {"_id": "$type", "count": {"$sum": 1}}}]))
        for x in data:
            h.write(str(x["_id"]) + "\t" + str(x["count"]) + "\n")
    if os.path.exists(stats_dir + "tax_tool.tsv"):
        os.remove(stats_dir + "tax_tool.tsv")
    shutil.move("tax_tool.tsv", stats_dir)

    with open("tax_barcode.tsv", "w") as h:
        h.write("#query\t#taxID\n")
        for x in db.barcodes.find({"tax": {"$exists": 1}}, {"processid": 1, "tax": 1}):
            h.write(x["processid"] + "\t" + str(int(x["tax"])) + "\n")
    execute('ktImportTaxonomy  -n %s -o  tax_barcode.html tax_barcode.tsv' % "Barcodes")
    if os.path.exists(stats_dir + "tax_barcode.html"):
        os.remove(stats_dir + "tax_barcode.html")
    shutil.move("tax_barcode.html", stats_dir)

    taxmap = {}

    with open("tax_genome.tsv", "w") as h:
        h.write("#query\t#taxID\n")
        for x in db.sequence_collection.find({"tax": {"$exists": 1}}):
            taxmap[x["name"]] = int(x["tax"]["tid"])
            h.write(x["name"] + "\t" + str(int(x["tax"]["tid"])) + "\n")
    execute('ktImportTaxonomy -n %s -o  tax_genome.html tax_genome.tsv' % "Genomas")
    if os.path.exists(stats_dir + "tax_genome.html"):
        os.remove(stats_dir + "tax_genome.html")
    shutil.move("tax_genome.html", stats_dir)

    with open("tax_seq.tsv", "w") as h:
        h.write("#query\t#taxID\n")
        for x in db.contig_collection.find({}, {"organism": 1, "name": 1}):
            if x["organism"] in taxmap:
                h.write(x["name"] + "\t" + str(taxmap[x["organism"]]) + "\n")
    execute('ktImportTaxonomy -n %s -o  tax_seq.html tax_seq.tsv' % "Seqs_Ensambl")
    if os.path.exists(stats_dir + "tax_seq.html"):
        os.remove(stats_dir + "tax_seq.html")
    shutil.move("tax_seq.html", stats_dir)

    with open("tax_prot.tsv", "w") as h:
        h.write("#query\t#taxID\n")
        for x in db.proteins.find({}, {"name": 1, "organism": 1}):
            if x["organism"] in taxmap:
                h.write(x["name"] + "\t" + str(taxmap[x["organism"]]) + "\n")
    execute('ktImportTaxonomy -n %s -o  tax_prot.html tax_prot.tsv' % "Proteinas")
    if os.path.exists(stats_dir + "tax_prot.html"):
        os.remove(stats_dir + "tax_prot.html")
    shutil.move("tax_prot.html", stats_dir)

    with open("tax_struct.tsv", "w") as h:
        h.write("#query\t#taxID\n")
        for x in dbpdb.structures.find({"tax": {"$exists": 1}}, no_cursor_timeout=True):
            h.write(x["name"] + "\t" + str(x["tax"]) + "\n")
    execute('ktImportTaxonomy  -n %s -o  tax_struct.html tax_struct.tsv' % "Estructuras")
    if os.path.exists(stats_dir + "tax_struct.html"):
        os.remove(stats_dir + "tax_struct.html")
    shutil.move("tax_struct.html", stats_dir)
    # http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=info&id= replace  http://localhost:8080/sndg/search/results?type=genome&query=&pageSize=50&start=0&taxonomia=