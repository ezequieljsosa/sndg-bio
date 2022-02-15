"""

https://github.com/tseemann/prokka
"""
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
        bpio.write(contigs,hw_gb,"gb")

    def proteins(self, h_or_str_gb, h_or_str_faa, otype="prot"):
        if isinstance(h_or_str_gb, str):
            h_gb = smart_parse(h_or_str_gb)
        else:
            h_gb = h_or_str_gb

        if isinstance(h_or_str_faa, str):
            h_faa = open(h_or_str_faa, "w")
        else:
            h_faa = h_or_str_faa

        try:
            for contig in h_gb:
                org = (" [" + contig.annotations["organism"].replace("[", "_").replace("]", "_") + "]"
                       ) if "organism" in contig.annotations else ""
                for feature in contig.features:
                    if feature.type in ["CDS", "RNA", "mat_peptide"]:
                        seq = None
                        description = feature.qualifiers["product"][0] if "product" else (
                            feature.qualifiers["note"][0] if "note" in feature.qualifiers else "")

                        locus_tag = feature.qualifiers["locus_tag"][0]
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
                            bpio.write(record, h_faa, "fasta")
                        else:
                            assert "pseudo" in feature.qualifiers


        finally:
            if isinstance(h_or_str_faa, str):
                h_faa.close()
            try:
                h_gb.close()
            except:
                pass


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
    cmd.add_argument('output_dir', nargs='?', default=sys.stdout)
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
            gbk_h = fileinput.input(args.input_bgk)
        utils.proteins(gbk_h, args.output_dir,otype=args.otype)
    if args.command == "fix":
        utils.fixgb(fileinput.input(args.input_bgk), args.output_gb)