import gzip
from SNDG import execute
import os
from tqdm import tqdm
import Bio.SeqIO as bpio
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from exceptions import ValueError

FNULL = open(os.devnull, 'w')

if __name__ == '__main__':
    workdir = "/data/bio/patho2/"
    workdir2 = "/data/bio/patho2_results/"

    with tqdm(os.listdir(workdir)) as pbar:
        for x in pbar:
            pbar.set_description(x)

            if not os.path.isdir(workdir + x):
                continue
            if os.path.exists(workdir2 + x.replace(" ", "_") + ".faa"):
                continue
            with open(workdir2 + x.replace(" ", "_") + ".faa", "w") as proteome_handler:
                with open(workdir2 + x.replace(" ", "_") + ".fna", "w") as gene_handler:
                    for g in tqdm(os.listdir(workdir + x)):
                        if g.startswith("GCF"):
                            with gzip.open(workdir + x + "/" + g) as assembly_handle:
                                try:
                                    for contig in bpio.parse(assembly_handle, "gb"):
                                        for f in contig.features:
                                            if (f.type == "CDS") and ("pseudo" not in f.qualifiers) and ("translation" in f.qualifiers) and ("locus_tag" in f.qualifiers):
                                                lt = f.qualifiers["locus_tag"][0]
                                                desc = "||".join([y + "=" + f.qualifiers[y][0]
                                                                  for y in ["product", "old_locus_tag", "protein_id"]
                                                                  if y in f.qualifiers])
                                                r = SeqRecord(id=lt, name=lt, description=desc,
                                                              seq=Seq(f.qualifiers["translation"][0]))
                                                bpio.write(r, proteome_handler, "fasta")

                                                r = SeqRecord(id=lt, name=lt, description=desc,
                                                              seq=f.extract(contig.seq))
                                                bpio.write(r, gene_handler, "fasta")
                                except ValueError:
                                    with open(workdir + "errors.txt", "a") as h:
                                        h.write(x + "\t" + g + "\n")
