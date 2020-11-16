import os
from tqdm import  tqdm

class GFFUtils:

    @staticmethod
    def gene_list(gff_file_path,gene_ftype="mRNA",info_gene_id="ID",tqdm_fn=tqdm):
        genes = {}
        with open(gff_file_path) as h:
            if tqdm_fn:
                h = tqdm_fn(h)
            for line in h:
                if line .startswith("#"):
                    continue

                vec = line.strip().split("\t")
                contig,_,ftype = vec[0:3]
                if ftype == gene_ftype:
                    info = {x.split("=")[0]:x.split("=")[1] for x in vec[-1].split(";")}
                    genes[ info[info_gene_id]] = 1

        return list(genes)




