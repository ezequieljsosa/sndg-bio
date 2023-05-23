#!/usr/bin/env python3

import random
import sys

import Bio.SeqIO as bpio
from Bio.SeqRecord import SeqRecord
import logging

_log = logging.getLogger(__name__)


class FastaSampler:

    def __init__(self, path_fasta, read_size=250, step=100, depth=100):
        self.path_fasta = path_fasta
        self.read_size = read_size
        self.step = step
        self.depth = depth
        self.seq_dict = None

    def parse(self):
        self.seq_dict = bpio.to_dict(bpio.parse(self.path_fasta, "fasta"))

    def total(self):
        total = 0
        for seq in self.seq_dict.values():
            total += int(len(seq) / self.step) * self.depth
        return total

    def create_pair_fastq(self):
        read_id = 0
        for c, seq in self.seq_dict.items():
            seq_len = len(seq)
            for curr_step in range(int(seq_len / self.step) - 1 ):
                pos = curr_step * self.step
                for i in range(self.depth):
                    read_id = read_id + 1
                    pos_noise = random.randint(-30, 30)
                    size_noise = random.randint(-10, 0)
                    pos_and_noise = pos + pos_noise
                    pos_and_noise = pos_and_noise if pos_and_noise > 0 else 0
                    subseq = seq.seq[pos_and_noise:pos_and_noise + self.read_size + size_noise]

                    assert subseq, [pos_and_noise, pos_and_noise + self.read_size + size_noise]

                    record = SeqRecord(id=f'seq_id_{str(read_id).zfill(20)}', description="",
                                       seq=subseq,
                                       letter_annotations={
                                           "phred_quality": [60 for _ in range(len(subseq))]})

                    subseq2 = seq.seq[pos_and_noise+ self.read_size:pos_and_noise + 2 * self.read_size + size_noise]

                    if len(subseq2):
                        assert subseq2, [pos_and_noise+ self.read_size,pos_and_noise + 2 * self.read_size + size_noise,subseq2]

                        record2 = SeqRecord(id=f'seq_id_{str(read_id).zfill(20)}', description="",
                                       seq=subseq2,
                                       letter_annotations={
                                           "phred_quality": [60 for _ in range(len(subseq2))]})
                        yield record,record2

    def create_single_fastq(self):
        read_id = 0
        for c, seq in self.seq_dict.items():
            seq_len = len(seq)
            for curr_step in range(int(seq_len / self.step)):
                pos = curr_step * self.step
                for i in range(self.depth):
                    read_id = read_id + 1
                    pos_noise = random.randint(-30, 30)
                    size_noise = random.randint(-10, 0)
                    pos_and_noise = pos + pos_noise
                    pos_and_noise = pos_and_noise if pos_and_noise > 0 else 0
                    subseq = seq.seq[pos_and_noise:pos_and_noise + self.read_size + size_noise]

                    assert subseq, [pos_and_noise, pos_and_noise + self.read_size + size_noise]

                    record = SeqRecord(id=f'seq_id_{str(read_id).zfill(20)}', description="",
                                       seq=subseq,
                                       letter_annotations={
                                           "phred_quality": [60 for _ in range(len(subseq))]})
                    yield record


if __name__ == "__main__":
    import os
    import argparse
    from tqdm import tqdm

    parser = argparse.ArgumentParser(description='Creates a FastQ file from a fasta file')

    parser.add_argument('--path_fasta', action='store', type=str,  help="fasta file")
    parser.add_argument('--read_size', action='store', type=int,
                        default=250)
    parser.add_argument('--step', action='store', type=int,
                        default=100, help="base pairs between the start position of each read batch")
    parser.add_argument('--depth', action='store', type=int,
                        default=50)

    parser.add_argument('-o', '--output', action='store', type=str, default=None)
    parser.add_argument('-o2', '--output2', action='store', type=str, default=None)

    parser.add_argument('-v', '--verbose', action="store_true")
    parser.add_argument('-s', '--silent', action="store_true")

    args = parser.parse_args()

    if not args.verbose:
        if os.environ.get('VERBOSE'):
            args.verbose = True

    if args.silent:
        _log.disabled = True

    assert os.path.exists(args.path_fasta), f"'{args.path_fasta}' does not exist"

    fs = FastaSampler(args.path_fasta, args.read_size, args.step, args.depth)
    fs.parse()

    h = sys.stdout
    if args.output:
        h = open(args.output,"w")
    if args.output2:
        with open(args.output2,"w") as h2:
            for r1,r2 in tqdm(fs.create_pair_fastq(), total=fs.total(), file=sys.stderr):
                bpio.write(r1, h, "fastq-illumina")
                bpio.write(r2, h2, "fastq-illumina")
    else:
        for r in tqdm(fs.create_single_fastq(), total=fs.total(), file=sys.stderr):
            bpio.write(r, h, "fastq-illumina")

    if args.output:
        h.close()
