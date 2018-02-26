"""

"""
from collections import defaultdict



import Bio.SeqIO as bpio
import Bio.SearchIO as bpsio
from BCBio import GFF


def identity(hsp):
    return 1.0 * hsp.ident_num / hsp.aln_span


def coverage(query_result, hsp):
    # return 1.0 * hsp.aln_span / hsp.query_span
    return 1.0 * hsp.query_span / query_result.seq_len


def hit_coverage(hit, hsp):
    # return 1.0 * hsp.aln_span /  hsp.hit_span
    return 1.0 * hsp.hit_span / hit.seq_len


blast_columns = ["query", "hit", "identity", "aln_len", "mismatches", "gap_openings", "qstart", "qend", "hstart",
                 "hend", "evalue"]


def read_blast_table(blast_table_file):
    import pandas as pd
    return pd.read_table(blast_table_file,
                  sep="\t", names=blast_columns, index_col=False)


def search_iterator(bps_iter):
    for query in bps_iter:
        for hit in query:
            for hsp in hit:
                yield (query, hit, hsp)


def add_blast_xml_props(sndg_iter):
    for query, hit, hsp in sndg_iter:
        hsp.identity = identity(hsp)
        hsp.coverage = coverage(query, hsp)
        hsp.hit_coverage = hit_coverage(hit, hsp)
        yield (query, hit, hsp)


def smart_parse(path):
    if path.endswith(".fasta"):
        return bpio.parse(path, "fasta")
    if path.endswith(".faa"):
        return bpio.parse(path, "fasta")
    if path.endswith(".fna"):
        return bpio.parse(path, "fasta")

    if path.endswith(".gb"):
        return bpio.parse(path, "gb")
    if path.endswith(".gbk"):
        return bpio.parse(path, "gb")
    if path.endswith(".genebank"):
        return bpio.parse(path, "gb")
    if path.endswith(".gbff"):
        return bpio.parse(path, "gb")

    if path.endswith(".embl"):
        return bpio.parse(path, "embl")

    if path.endswith(".gbff"):
        with open(path) as h:
            return GFF(h)

    if path.endswith(".fq"):
        return bpio.parse(path, "fastq")
    if path.endswith(".fastq"):
        return bpio.parse(path, "fastq")

    if path.endswith(".hmm"):
        return search_iterator(bpsio.parse(path, "hmmer3-text"))

    if path.endswith(".xml"):
        with open(path) as h:
            h.readline()
            l = h.readline()
            if "BlastOutput" in l:
                return add_blast_xml_props(search_iterator(bpsio.parse(path, "blast-xml")))
            if "<uniprot" in l:
                return bpio.parse(path, "uniprot-xml")

    raise Exception("invalid format")


def blast_parse(path, fn_hsp_filter):
    for query, hit, hsp in add_blast_xml_props(search_iterator(bpio.parse(path, "uniprot-xml"))):
        if fn_hsp_filter(query, hit, hsp):
            yield (query, hit, hsp,)


def blast_partition(path, partition_dict):
    partitions = defaultdict(lambda: [], {"no_hit": [], "else": []})
    for query in bpsio.parse(path, "blast-xml"):
        hits = list(query)
        if not hits:
            partitions["no_hit"].append(query.id)
        for hit in hits:
            for hsp in hit:
                hsp.identity = identity(hsp)
                hsp.coverage = coverage(query, hsp)
                hsp.hit_coverage = hit_coverage(hit, hsp)

                added = False
                for k, fn_filter in partition_dict.items():
                    if fn_filter(query, hit, hsp):
                        partitions[k].append((query, hit, hsp,))
                        added = True
                        break
                if not added:
                    partitions["else"].append((query, hit, hsp,))

    return partitions


if __name__ == '__main__':
    # next(smart_parse("/data/pext/good_hits.xml"))
    blast_partition("/tmp/GCF_001624625.1/pdb_blast.xml",
                    {"cristal": lambda q, h, hsp: (1.0 * hsp.ident_num / hsp.aln_span) > 0.9})
