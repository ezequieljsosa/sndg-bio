

class AnnMerge:

    def merge(self,source_contigs:dict,target_contigs:list,hw,locus_tag:str):
        """

        :param source_contigs: dict of contigs with features to add
        :param target_contigs: list of contigs to add the source features
        :return:
        """
        for target_contig in target_contigs:
            if target_contig.id in source_contigs:
                for source_feature in source_contigs[target_contig.id]:

