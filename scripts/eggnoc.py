import logging
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

from SNDG import init_log,execute
import  os

init_log()

_log = logging.getLogger("eggnoc")

from multiprocessing import Process

if __name__ == "__main__":
    parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-dir", "--faas_dir", required=True)
    parser.add_argument("-tax", "--tax_map", required=True)
    parser.add_argument("-o", "--out_dir", default="./")
    parser.add_argument( "--cpu", type=int, default=1)


    args = parser.parse_args()

    args.out_dir = os.path.abspath(args.out_dir)

    tax_map = {l.strip().split("\t")[0]:l.strip().split("\t")[1] for l in open(args.tax_map)}
    if "bact" in tax_map.values():
        def server():
            execute("python emapper.py -d bact --cpu %i --servermode" % args.cpu)
        p = Process(target=server)
        p.start()

    for accession,tax in tax_map.items():
        if tax == "bact":
            execute("python emapper.py -d bact:localhost:51600 -i %s.faa -o %s" %
                    (accession ,accession) )
        else:
            execute("python emapper.py -i %s.faa --output %s -d euk" %
                    (args.out_dir + accession , args.out_dir + accession) )


