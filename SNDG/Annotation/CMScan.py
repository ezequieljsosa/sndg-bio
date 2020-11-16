import os
from BCBio import GFF
from Bio.SeqFeature import SeqFeature, FeatureLocation


class CMScan():
    """
    sed  's/  */  /g' genome.tblout > ./genome.tblout2
    grep -v rRNA genome.tblout2 | grep -v "^#" | grep -v tRNA > other_rnas.tbl
     cmscan --cut_ga --rfam --nohmmonly --tblout genome.tblout --fmt 2 --clanin /data/Rfam.clanin /data/Rfam.cm ref.fa.masked > genome.cmscan
    # test --cut_tc
    """

    def __init__(self, file_path):
        self.file_path = file_path
        self.records = []

    def load(self):
        fields = ["hit", "accession", "contig", "start", "end", "strand", "desc"]
        with open(self.file_path) as h:
            for line in h:
                r = line.strip().split("\t")
                data = [r[1], r[2], r[3], r[9], r[10], r[11], " ".join(r[26:]).replace(";",",") ]
                record = {fields[i]: data[i] for i in range(len(fields))}
                self.records.append(record)

    def add2gff(self, gff_in, gff_out, locus_tag="otherRNA_"):
        rs = []
        idx = 1
        with open(gff_in) as h:
            for seqrecord in GFF.parse(h):
                rs.append(seqrecord)
                records = [x for x in self.records if x["contig"] == seqrecord.id]
                for record in records:
                    rid = locus_tag + str(idx)
                    idx += 1
                    loc = FeatureLocation(start=int(record["start"]), end=int(record["end"]),
                                          strand=1 if record["strand"] == "+" else -1)
                    sf = SeqFeature(location=loc,
                                    type=record["hit"], qualifiers={
                            "Parent": rid, "RFAM": record["accession"],
                            "product": record["desc"],
                            "ID": rid + "-T1", "locus_tag": rid + "-T1"})
                    f = SeqFeature(location=loc,
                                   type="gene", qualifiers={"ID": rid, "locus_tag": rid}
                                   )
                    f.sub_features=[sf]
                    seqrecord.features.append(f)
        with open(gff_out, "w") as h:
            GFF.write(rs, h)
