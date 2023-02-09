"""

https://github.com/tseemann/prokka
"""
from Bio.SeqFeature import SeqFeature, FeatureLocation

from SNDG import docker_wrap_command, DOCKER_MAPPINGS, execute
import fileinput
from SNDG.Sequence import smart_parse
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import Bio.SeqIO as bpio
import sys
from tqdm import tqdm


class GenebankUtils:
    """

    """

    def __init__(self):
        pass

    def fixgb(self, h_gb, hw_gb):
        contigs = list(bpio.parse(h_gb, "gb"))
        # for contig in contigs:
        #     features = []
        #     for f in contig:
        #         if f.type in ["CDS","mRNA"]

        for contig in contigs:
            for f in contig.features:
                if f.type == "CDS":
                    if "translation" not in f.qualifiers:
                        seq = str(f.extract(contig).seq.translate())
                        if seq.count("*") > 1:
                            f.qualifiers["pseudo"] = ["true"]
                        else:
                            f.qualifiers["translation"] = [seq]
        bpio.write(contigs, hw_gb, "gb")

    def lt(self, feature):
        return feature.qualifiers.get("locus_tag", [""])[0]

    def region_from_lts(self, lts, contig, region_name, rtype="REGION", strand_lts=None):

        region = []
        strand = None

        for f in contig.features:
            if self.lt(f) in lts:
                region.append(f)
                if strand_lts and (self.lt(f) in strand_lts):
                    strand = f.location.strand

        locations = [x.location.start for x in region] + [x.location.end for x in region]
        if locations:
            return SeqFeature(type=rtype, location=FeatureLocation(
                start=min(locations), end=max(locations), strand=strand), qualifiers={"name": [region_name]})
        return None

    def proteins(self, sequences, otype="prot"):

        for contig in sequences:
            org = (" [" + contig.annotations["organism"].replace("[", "_").replace("]", "_") + "]"
                   ) if "organism" in contig.annotations else ""
            for feature in contig.features:
                if feature.type in ["CDS", "RNA", "mat_peptide"]:
                    seq = None
                    location = f'{contig.id}:{feature.location.start}-{feature.location.end}'
                    description = location + " " + feature.qualifiers["product"][0] if "product" else (
                        feature.qualifiers["note"][0] if "note" in feature.qualifiers else "")

                    locus_tag = feature.qualifiers.get("locus_tag", [location + "_" + feature.type])[0]
                    gene = feature.qualifiers["gene"][0] if "gene" in feature.qualifiers else ""

                    if feature.type == "mat_peptide":
                        gene = gene + "_" + feature.qualifiers["product"][0]
                        locus_tag = locus_tag + "_" + feature.qualifiers["product"][0].replace(" ", "_")
                    if otype == "prot":
                        if feature.type == "CDS" and "pseudo" not in feature.qualifiers:
                            seq = Seq(feature.qualifiers["translation"][0])
                        elif feature.type == "mat_peptide":
                            seq = feature.extract(contig.seq).translate()
                    else:
                        seq = feature.extract(contig.seq)
                    if seq:
                        record = SeqRecord(id=locus_tag, name=gene, description=description + org, seq=seq)
                        yield record
                    else:
                        assert "pseudo" in feature.qualifiers

        # finally:
        #     if isinstance(h_or_str_faa, str):
        #         h_faa.close()
        #     try:
        #         h_gb.close()
        #     except:
        #         pass


if __name__ == '__main__':
    import argparse
    import os
    from SNDG.Sequence import smart_parse

    parser = argparse.ArgumentParser(description='Utils over genebank file')

    subparsers = parser.add_subparsers(help='commands', description='valid subcommands', dest='command', required=True)
    cmd = subparsers.add_parser('fix', help='fix genebank file to be imported')
    cmd.add_argument('input_bgk')
    cmd.add_argument('output_gb', nargs='?', default=sys.stdout)

    cmd = subparsers.add_parser('genes', help='extracts the list of gene products from ')
    cmd.add_argument('input_bgk')
    cmd.add_argument('output_fna', nargs='?', default=sys.stdout)
    cmd.add_argument('-otype', help="output type", choices=["prot", "nucl"], default="prot")

    args = parser.parse_args()
    utils = GenebankUtils()

    if args.command == "genes":
        if isinstance(args.input_bgk, str):
            if not os.path.exists(args.input_bgk):
                sys.stderr.write(f'{args.input_bgk} not found')
                sys.exit(1)
        if args.input_bgk != "-":
            gbk_h = smart_parse(args.input_bgk)
        else:
            gbk_h = bpio.parse(fileinput.input(args.input_bgk), "gb")

        if isinstance(args.output_fna, str):
            h_faa = open(args.output_fna, "w")
        else:
            h_faa = args.output_fna
        try:
            for seq in utils.proteins(gbk_h, otype=args.otype):
                bpio.write(seq, h_faa, "fasta")
        finally:
            try:
                h_faa.close()
            except:
                pass

    if args.command == "fix":
        utils.fixgb(fileinput.input(args.input_bgk), args.output_gb)
