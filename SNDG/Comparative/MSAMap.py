import fileinput
from collections import defaultdict

import Bio.SeqIO as bpio
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class MSAMap():

    @staticmethod
    def from_msa(msa_file_or_path):
        seqs = {x.id: x for x in bpio.parse(fileinput.input(msa_file_or_path), "fasta")}
        msa = MSAMap(seqs)
        msa.init()
        return msa

    def __init__(self, seqs: dict, gap_code="-"):
        self.seqs = seqs
        self.pos_msa_seq_map = {x: {} for x in self.seqs.keys()}
        self.pos_seq_msa_map = {x: {} for x in self.seqs.keys()}
        self.gap_code = gap_code

    def aln_len(self):
        return len(self.seqs[list(self.seqs)[0]])

    def samples(self):
        return list(self.seqs.keys())

    def init(self):
        curr_pos = {x: 0 for x in self.samples()}
        for msa_pos in range(self.aln_len()):
            for sample in self.samples():
                if self.seqs[sample][msa_pos] != self.gap_code:
                    self.pos_msa_seq_map[sample][msa_pos] = curr_pos[sample]
                    self.pos_seq_msa_map[sample][curr_pos[sample]] = msa_pos
                    curr_pos[sample] += 1

    def variants(self, refseq):
        variants = defaultdict(lambda: defaultdict(list))
        for msa_pos in range(self.aln_len()):
            if msa_pos in self.pos_msa_seq_map[refseq]:
                ref_pos = self.pos_msa_seq_map[refseq][msa_pos]
                ref_data = self.seqs[refseq][msa_pos]
                assert ref_data != "-"
                variant_id = f"{ref_data}_{ref_pos}"
                for sample in self.samples():
                    data = self.seqs[sample][msa_pos]
                    if ref_data == self.gap_code:
                        variants[variant_id]["del" + data].append(sample)
                    elif ref_data != data:
                        variants[variant_id][data].append(sample)
            else:
                for sample in self.samples():
                    data = self.seqs[sample][msa_pos]
                    if data != self.gap_code:
                        variants[variant_id]["ins" + data].append(sample)

        return {k: {k1: v1 for k1, v1 in v.items()} for k, v in variants.items()}

    def pos_from_seq(self, seq_name_in, pos_in, seq_name_out):
        if pos_in in self.pos_seq_msa_map[seq_name_in]:
            msa_pos = self.pos_seq_msa_map[seq_name_in][pos_in]
            if msa_pos in self.pos_msa_seq_map[seq_name_out]:
                return self.pos_msa_seq_map[seq_name_out][msa_pos]
            else:
                raise ValueError(f"MSAPos {msa_pos} from {seq_name_out} not found")
        else:
            raise ValueError(f"Pos {pos_in} from {seq_name_in} not found")

    def exists_pos(self, seq_name_in, pos_in, seq_name_out):
        return (pos_in in self.pos_seq_msa_map[seq_name_in]) and (
                self.pos_seq_msa_map[seq_name_in][pos_in] in self.pos_msa_seq_map[seq_name_out])

    def subseq(self, seq_name_in, pos_in_start, pos_in_end, seq_name_out):
        seq = ""
        start = -1

        for pos_in in range(pos_in_start, pos_in_end):
            try:
                pos_out = self.pos_from_seq(seq_name_in, pos_in, seq_name_out)
                seq += self.seqs[seq_name_out][pos_out]
                if start == -1:
                    start = pos_out
            except ValueError:
                continue
        end = pos_out

        return SeqRecord(id=f'{seq_name_out}_{start}_{end}', name="", description="", seq=Seq(seq))

    def vcf(self, ref_sequence, stdout=sys.stdout):
        assert ref_sequence in self.samples(), f'{ref_sequence} does not exists'
        samples = ref_sequence + sorted(set(self.samples()) - set([ref_sequence]))

        header = ("\t".join("#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT".split()) +
                  "\t" + "\t".join(samples))
        stdout.write(f"""##fileformat=VCFv4.2
    ##ALT=<ID=*,Description="Represents allele(s) other than observed.">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    {header}\n""")

        first_pos = 0
        pos_ant = -1
        acc_ref = None
        acc_alt = None
        for x, v in msa.variants(args.ref_sequence).items():
            ref, pos = x.split("_")
            if not first_pos:
                first_pos = pos
            alts = {alt: i + 1 for i, alt in enumerate(v)}
            sample_dict = {}
            for alt, alt_samples in v.items():
                for sample in alt_samples:
                    sample_dict[sample] = alts[alt]
            sample_gts = "0\t" + "\t".join(
                [str(sample_dict[sample]) for sample in samples if sample != args.ref_sequence])
            alts_str = ','.join([str(x) for x in alts.keys()]).upper()

            if (int(pos) - 1) == int(pos_ant):
                if acc_ref:
                    acc_ref += ref if ref.replace("-", "") else ""
                    acc_alt += alts_str if alts_str.replace("-", "") else ""
                else:
                    acc_ref = ref
                    acc_alt = alts_str
            else:
                if acc_ref:
                    if acc_alt == msa.gap_code:
                        # deletion
                        ant_ref = str(
                            msa.subseq(args.ref_sequence, int(first_pos), int(first_pos) + 1, args.ref_sequence).seq)
                        acc_ref = (ant_ref + acc_ref).upper()
                        acc_alt = ant_ref.upper()
                        first_pos = str(int(first_pos) - 1)
                    elif acc_ref == msa.gap_code:
                        # insertion
                        ant_ref = str(
                            msa.subseq(args.ref_sequence, int(first_pos), int(first_pos) + 1, args.ref_sequence).seq)
                        acc_ref = ant_ref.upper()
                        acc_alt = (ant_ref + acc_alt).upper()
                        first_pos = str(int(first_pos) - 1)
                    line = f"{args.ref_sequence}\t{int(first_pos) + 1}\t.\t{acc_ref.upper()}\t{acc_alt}\t.\t.\t.\tGT\t{sample_gts}"
                    stdout.write(line + "\n")
                    acc_ref = None
                    acc_alt = None
                    first_pos = pos
                acc_ref = ref
                acc_alt = alts_str
            pos_ant = pos
        if acc_ref:
            line = f"{args.ref_sequence}\t{int(first_pos) + 1}\t.\t{acc_ref.upper()}\t{acc_alt}\t.\t.\t.\tGT\t{sample_gts}"
            stdout.write(line + "\n")


if __name__ == '__main__':
    import argparse
    import os

    parser = argparse.ArgumentParser(description='Mapping to variant calls pipeline.')
    required = parser.add_argument_group('required arguments')
    required.add_argument('-i', '--input_msa', action='store', default="-")
    required.add_argument('-r', '--ref_sequence', action='store', required=True)
    required.add_argument('-o', '--output', default="-")

    args = parser.parse_args()

    if args.input_msa != "-" and not os.path.exists(args.input_msa):
        raise FileNotFoundError(f"{args.input_msa} does not exists")

    seqs = {x.id: x for x in bpio.parse(fileinput.input(args.input_msa), "fasta")}
    msa = MSAMap(seqs)
    msa.init()

    if args.output == "-":
        msa.vcf(args.ref_sequence, stdout=sys.stdout)
    else:
        with open(args.output, "w") as h:
            msa.vcf(args.ref_sequence, stdout=h)
