import vcf
from vcf.model import _Call,make_calldata_tuple

vcf_reader = vcf.Reader(open('./arch_final.vcf', 'r'))
h = open('./pepe.vcf', 'w')
vcf_writer = vcf.Writer(h, vcf_reader)

for record in vcf_reader:
            if record.ALT and record.ALT[0]:
                dcc = make_calldata_tuple(["GT"])
                ngt = len(record.alleles)
                dc = dcc(str(ngt))
                record.alleles.append("X")
                record.ALT.append("X")
                c = _Call(record,"ANNPOS",dc)
                record.samples.append(c)
                vcf_writer.write_record(record)
                break
vcf_writer.close()
