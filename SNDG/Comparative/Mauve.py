class Mauve():

    @classmethod
    def parse_orthologs(cls, file_path):
        orthologs = []
        with open(file_path) as h:
            for line in h:
                if line.strip():
                    ortho_dict = {ortho.split(":")[0]: ortho.split(":")[1] for ortho in line.strip().split("\t")}
                    orthologs.append(ortho_dict)
        return orthologs

    @classmethod
    def count_orthologs(cls, parsed_orthologs, ref_num="0"):
        """

        :param parsed_orthologs: result of Mauve.parse_orthologs
        :return:
        """
        count = {}
        for ortho in parsed_orthologs:
            if ref_num in ortho:
                count[ortho[ref_num]] = len(ortho)
        return count


if __name__ == '__main__':
    from SNDG import init_log
    from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB

    init_log()
    mdb = BioMongoDB("tdr",port=27018)
    datafile = "/data/organismos/SaureusN315/annotation/conservation/target_props.tsv"
    parsed_orthologs = Mauve.parse_orthologs(
        "/data/organismos/SaureusN315/annotation/conservation/ortologos_staphylo.csv")
    count = Mauve.count_orthologs(parsed_orthologs, "0")
    with open(datafile, "w") as h:
        h.write("id\tconserved_count\tconserved_percent\n")
        max_count = max(count.values())
        for gene, count in count.items():
            h.write(gene + "\t" + str(count) + "\t" + ("%0.2f" % (count * 1.0 / max_count)) + "\n")

    mdb.load_metadata("SaureusN315",datafile)
