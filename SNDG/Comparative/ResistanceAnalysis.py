import pandas as pd
from tqdm import tqdm
import subprocess
from SNDG.Comparative.VcfSnpeffIO import VcfSnpeffIO


class ResistanceAnalysis:
    def __init__(self, path_db="/data/databases/tbprofiler.tbl"):
        self.resist = pd.read_csv(path_db)
        self.resist["AApos"] = [int(x) if x != "-" else "" for x in self.resist.AApos]
        self.resist["LocusTag"] = [x.split("_")[0] for x in self.resist.GeneID]
        self.resist["NucleotidePosH37"] = [int(x.split("/")[0]) if x != "-" else "" for x in
                                           self.resist.NucleotidePosH37]

        self.result = []

    # TODO filtro marcadores filogeneticos
    # TODO guardar tabla final
    # TODO Probar variables independientes (CAC	GAG) o (GG	CA)
    # TODO ver variable card +nt420:GG

    def process_vcf(self, vcf_path):
        variants = VcfSnpeffIO.parse(vcf_path)
        total = int(subprocess.check_output('grep -vc "^#" ' + vcf_path, shell=True))
        reported_genes = set(self.resist["LocusTag"])

        out_genes = []

        for variant, effects in tqdm(variants, total=total):

            main_eff = effects[0]

            # ver si esta en un gen reportado
            if main_eff.geneid in reported_genes:

                reported_gene_variants = self.resist[self.resist["GeneID"] == main_eff.geneid]
                out_genes.append(main_eff.geneid)
                # ver si la posicion esta reportada

                if main_eff.hgvs_p:
                    self.__process_protein_variant(main_eff, reported_gene_variants, variant)
                else:
                    self.__process_nucleotide_variant(main_eff, reported_gene_variants, variant)

        return self.result

    def __process_nucleotide_variant(self, main_eff, reported_gene_variants, variant):
        reported_positions = reported_gene_variants[
            reported_gene_variants["NucleotidePosH37"] == variant.POS]
        if not reported_positions.empty:

            reported = reported_positions[(reported_positions.REF == variant.REF) &
                                          (reported_positions.ALT.isin([str(x) for x in variant.ALT]))]
            if reported.empty:
                drug = set(reported_positions["Drug"])
            else:
                drug = set(reported["Drug"])

            for sample in variant.samples:
                if sample.called:
                    ad = sample.data.AD if hasattr(sample.data, "AD") else []
                    resist_record = {
                        "sample": sample.sample,
                        "gene": main_eff.geneid,
                        "depth": ad,
                        "drug": drug,
                        "contig": variant.CHROM,
                        "pos": variant.POS,
                        "c.hgvs": main_eff.hgvs_c,
                        "p.hgvs": None,
                        "reported": not reported.empty
                    }
                    self.result.append(resist_record)

    def __process_protein_variant(self, main_eff, reported_gene_variants, variant):
        # print (main_eff.aa_pos,list(reported_gene_variants["AApos"]))

        reported_positions = reported_gene_variants[reported_gene_variants["AApos"] == main_eff.aa_pos]
        if not reported_positions.empty:
            reported = reported_positions[(reported_positions.AAref == main_eff.aa_ref) &
                                          (reported_positions.AAalt == main_eff.aa_alt)]
            if reported.empty:
                drug = set(reported_positions["Drug"])
            else:
                drug = set(reported["Drug"])

            for sample in variant.samples:
                if sample.called:
                    ad = sample.data.AD if hasattr(sample.data, "AD") else []
                    resist_record = {
                        "sample": sample.sample,
                        "gene": main_eff.geneid,
                        "depth": ad,
                        "drug": drug,
                        "contig": variant.CHROM,
                        "pos": variant.POS,
                        "c.hgvs": main_eff.hgvs_c,
                        "p.hgvs": main_eff.hgvs_p,
                        "reported": not reported.empty
                    }
                    self.result.append(resist_record)


if __name__ == '__main__':
    import argparse
    from SNDG import init_log

    init_log()

    parser = argparse.ArgumentParser(description='Annotates the reported variants')
    required = parser.add_argument_group('required arguments')
    required.add_argument('-vcf', required=True)
    required.add_argument('-db', default="/data/databases/tbprofiler.tbl")

    args = parser.parse_args()

    ra = ResistanceAnalysis(args.db)
    data = ra.process_vcf(args.vcf)
    print len(data)
