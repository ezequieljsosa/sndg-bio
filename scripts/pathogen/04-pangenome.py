
import logging
import gzip

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
from glob import glob
from SNDG import init_log
import Bio.SeqIO as bpio

init_log(rootloglevel=logging.INFO)
_log = logging.getLogger(__name__)

if __name__ == "__main__":
    workdir = "/data/bartonella/ncbi-genomes-2018-06-25/"
    try:
        proth = open(workdir + "proteins.fasta","w")
        genesh = open(workdir + "genes.fasta","w")
        # pseudo = open(workdir + "proteins.fasta","w")
        for genebank_gz in tqdm(glob(workdir + "*.gz")):
            assembly = genebank_gz.split(workdir)[1].split("_genomic.gbff.gz")[0]
            for contig in tqdm( bpio.parse(gzip.open( genebank_gz),"gb") ):
                for f in contig.features:
                    if f.type == "CDS":

                        if "translation" in f.qualifiers:
                            locus_tag = f.qualifiers["locus_tag"][0]
                            desc = f.qualifiers["protein_id"][0] + " " + assembly
                            seq = SeqRecord(id=locus_tag, seq=Seq(f.qualifiers["translation"][0]),
                                            description= desc,name=locus_tag)
                            bpio.write(seq,proth,"fasta")
                            seq = f.extract(contig)
                            seq.id = locus_tag
                            seq.name = locus_tag
                            seq.descripion = desc
                            bpio.write(seq,genesh,"fasta")


    finally:
        proth.close()
        genesh.close()

    #cd-hit -c 0.9 -i proteins.fasta -o proteins_90.fasta  -g 1 -aS 0.8 -p 1
    #cd-hit -c 0.5 -i proteins.fasta -o proteins_50.fasta  -g 1 -aS 0.8 -p 1
    #cd-hit -c 0.9 -i genes.fasta -o genes_90.fasta  -g 1 -aS 0.8 -p 1
    #cd-hit -c 0.5 -i genes.fasta -o genes_50.fasta  -g 1 -aS 0.8 -p 1



