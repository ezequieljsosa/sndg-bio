import Bio.SeqIO as bpio
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='transforms Modpipe pdball.pir to a fasta format')
    parser.add_argument('-i', '--input', default="pdball.pir")
    parser.add_argument('-o', '--output', default="modeller.fasta")

    args = parser.parse_args()
    data = open(args.input).read()

    with open(args.output,"w") as h:
        for record in data.split("C; Produced by MODELLER"):
            record = [line.strip() for line in record.split("\n") if line.strip()]
            if record:
                seqrecord = SeqRecord(id=record[0][1:].strip(),
                                      name="",description=record[1].strip(),
                                      seq=Seq( "".join(record[2:])  ))
                bpio.write(seqrecord,h,"fasta")



