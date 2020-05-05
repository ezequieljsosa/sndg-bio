from tqdm import tqdm


class InterPro():

    @staticmethod
    def gff_get_go(gff_file_path, info_field="Ontology_term", tqdm_fn=tqdm):

        with open(gff_file_path) as h:
            if tqdm_fn:
                h = tqdm_fn(h)
            last_gene = ""
            gos = []
            for line in h:
                if line.startswith("#"):
                    continue
                vec = line.strip().split("\t")
                gene, _, ftype = vec[0:3]

                if (not last_gene) or (last_gene != gene):
                    yield (last_gene, gos)
                    gos = []
                    last_gene = gene


                info = {x.split("=")[0] if "=" in x else x: x.split("=")[1] if "=" in x else x for x in
                        vec[-1].split(";")}

                if info_field in info:
                    gos += info[info_field].replace('"', '').split(",")

        yield (last_gene, gos)
