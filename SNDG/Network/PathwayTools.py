"""
https://ask.pathwaytools.com/question/17972/access-pathologic-via-api/
So I defined organism-params.dat and genetic-elements.dat in a folder called dir_with_files, included the annotation file test.gbk and started pathwatools: pwt -patho dir_with_files

organism-params.dat
ID      test
NAME    Pseudomonas lurida
STORAGE File
NCBI-TAXON-ID   244566
DOMAIN  TAX-2

genetic-elements.dat
ID      TEST-CHROM-1
TYPE    :CHRSM
CIRCULAR?       Y
ANNOT-FILE      test.gbk

export PROXY=proxy.fcen.uba.ar:8080
pathway-tools -no-cel-overview -no-web-cel-overview  -patho /data/organismos/LHelv/annotation/pathways
"""
import logging
from tqdm import tqdm
from goatools.obo_parser import GODag

import Bio.SeqIO as bpio
from Bio import Alphabet
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from SNDG import execute, mkdir

_log = logging.getLogger(__name__)
go_dag = GODag("/data/databases/go/go.obo", load_obsolete=True)


def process_go_db(db, go_code):
    go = go_dag[go_code]
    if go.namespace == db and go.level > 0:
        return "|".join([go.name.strip().replace("\n", ""), go_code.split(":")[1], "1", "ISS"])
    else:
        return None


def process_go_cc(go_code):
    """
    GO term name, GO term ID, citation PubMed ID, and evidence code, separated by vertical bars. Example:
    antimicrobial humoral response|0019730|16163390|IMP
    http://geneontology.org/page/guide-go-evidence-codes -->
    Inferred from Sequence or structural Similarity (ISS) = default
    """
    return process_go_db("cellular_component", go_code)


def process_go_mf(go_code):
    return process_go_db("molecular_function", go_code)


def process_go_bp(go_code):
    return process_go_db("biological_process", go_code)


class PathwayTools:
    default_mappings = map_data = {
        "gene.mRNA": {
            "type": "CDS",
            "qualifiers": {
                "gene_symbol": "gene",
                "locus_tag": "locus_tag",
                "description": "product",
                "Note": "note",

                "EC": ("EC_number", lambda x: x.split(":")[1]),
                "x": "product_comment",
                "x": "gene_comment",
                "x": "pseudo",
                "x": "alt_name",
                "x": "db_xref",
                "GO": {"go_component": process_go_cc,
                       "go_function": process_go_mf,
                       "go_process": process_go_bp}
            }}
    }

    def __init__(self, workdir, assembly_level="CHRSM"):
        """
        :param workdir:
        :param assembly_level:  "CHRSM", "PLASMID", "MT" (mitochondrial chromosome),"PT" (chloroplast chromosome), or "CONTIG"
        """
        self.workdir = workdir
        self.assembly_level = assembly_level
        self.gb_file = workdir + "/pwtools.gb"

    def create_genetic_elements(self):
        """
        ;;==========================================================================
;; Sample genetic-elements.dat file, created 2/13/2001 by Suzanne Paley.
;; PathoLogic v. 1.8
;;
;; This file contains a set of records, one for each genetic element
;; (chromosome, plasmid, etc.) or contig in the organism.  Each genetic
;; element can be described either with a single annotation file or as a set
;; of contigs (but not both).  Contigs are each described with a single
;; annotation file.  Records are ended by a single line containing only the
;; string "//".  Each record contains a set of attribute-value pairs, one
;; per line, separated by the tab character.  Lines beginning with ; are
;; comment lines.  Valid attributes are:
;;
;; ID -- Required.  Unique identifier for the genetic element.  The identifer
;;       must be unique with respect to all other identifiers used in the
;;       PGDB, such as gene and protein identifiers.  The ID should start
;;       with an alphabetic character, and contain no spaces or special
;;       characters.  Case insensitive.
;; NAME -- Optional.  Descriptive name for the genetic element
;; TYPE -- Optional.  One of :CHRSM, :PLASMID, :MT (mitochondrial chromosome),
;;         :PT (chloroplast chromosome), or :CONTIG.  Defaults to :CHRSM.
;; CIRCULAR? -- Required (except for contigs).  Is the genetic element
;;              circular?  Valid options are Y or N.
;; CODON-TABLE -- Optional, defaults to the codon table specified in the
;;                organism.dat file, or to the standard codon-table if not
;;                supplied.  This should be a number between 1 and 15.
;; ANNOT-FILE -- Required, unless the CONTIG attribute is supplied.  Filename
;;               for the annotations file for this genetic element.  The file
;;               must be in either pathologic format (.pf) or genbank
;;               format (.gbk).  If the file is in the orgcyc/version/input/
;;               directory, then only the filename is required, otherwise a
;;               full pathname must be supplied.
;; SEQ-FILE -- Optional.  Filename for the sequence file for this
;;             genetic element.  The file must be in FASTA format
;;             (.fna or .fsa).  If the file is in the orgcyc/version/input/
;;             directory, then only the filename is required, otherwise a full
;;             pathname must be supplied.
;; CONTIG -- Optional (for all types other than :CONTIG.  This field is not
;;           permitted in records of type :CONTIG.)  ID for a contig that
;;           forms part of this genetic element.  A genetic element may have
;;           multiple CONTIG lines.  This field should not be supplied if the
;;           ANNOT-FILE attribute is supplied.
;;
;; Note that the "//" is required after the final genetic element, also.
;;===========================================================================

ID	TEST-CHROM-1
NAME	Chromosome 1
TYPE	:CHRSM
CIRCULAR?	N
ANNOT-FILE	chrom1.pf
SEQ-FILE	chrom1.fsa
//
ID	TEST-CHROM-2
NAME	Chromosome 2
CIRCULAR?	N
ANNOT-FILE	/mydata/chrom2.gbk
SEQ-FILE	/mydata/chrom2.fna
//
ID	TEST-CHROM-3
NAME	Chromosome 3
TYPE	:CHRSM
CIRCULAR?	N
CONTIG	CONTIG-1
CONTIG	CONTIG-2
//
ID	TEST-MIT-CHROM
NAME	Mitochondrial chromosome
TYPE	:MT
CIRCULAR?	Y
CODON-TABLE	2
ANNOT-FILE	mit-chrom.pf
//
ID	CONTIG-1
NAME	Contig 1 of Chromosome 3
TYPE	:CONTIG
ANNOT-FILE	chrom3-contig1.pf
SEQ-FILE	chrom3-contig1.fsa
//
ID	CONTIG-2
NAME	Contig 2 of Chromosome 3
TYPE	:CONTIG
ANNOT-FILE	chrom3-contig2.pf
SEQ-FILE	chrom3-contig2.fsa
//
        """
        mkdir(self.workdir + "/")
        with open(self.workdir + "/genetic-elements.dat", "w") as genetic_elements:
            for seq_record in bpio.parse(self.gb_file, "genbank"):
                new_gb = open(self.workdir + "/%s.gbk" % seq_record.id, "w")
                bpio.write(seq_record, new_gb, "genbank")

                genetic_elements.write("ID\t" + seq_record.id + "\n")
                genetic_elements.write("TYPE\t:" + self.assembly_level + "\n")
                genetic_elements.write("CIRCULAR?\tN" + "\n")
                genetic_elements.write("ANNOT-FILE\t" + seq_record.id + ".gbk" + "\n")
                genetic_elements.write("//" + "\n")

    def create_organism_params(self, name, organism, tax, domain):
        """
        :param name: Name of the collection in pwtools
        :param organism: Description of the organism
        :param tax: taxid
        :param domain: "TAX-2" (Bacteria), "TAX-2759" (Eukaryota), and "TAX-2157" (Archaea).
        :return:
        """
        template = """ID\t{name}
NAME\t{organism}
STORAGE\tFile
NCBI-TAXON-ID\t{tax}
DOMAIN\t{domain}"""
        with open(self.workdir + "organism-params.dat", "w")  as h:
            h.write(template.format(name=name, organism=organism, tax=tax, domain=domain))

    def create_genebank(self, contig_iterator, mappings=default_mappings):
        def dbfeature2seqfeature(org_feature):
            seqf = SeqFeature(
                FeatureLocation(org_feature.location.start, org_feature.location.end, org_feature.location.strand),
                org_feature.type, id=org_feature.id)

            self.map_atributes(mappings, org_feature, seqf)

            return seqf

        def process_contig(contig):
            record = SeqRecord(id=contig.id, seq=Seq(str(contig.seq), alphabet=Alphabet.DNAAlphabet()))
            for f in contig.features:

                if f.type.lower() == "gene":
                    seqfeature = dbfeature2seqfeature(f)
                elif f.type.lower() in ["rrna", "trna", "ncrna"]:
                    seqfeature = SeqFeature(f.location,
                                            f.type.replace("ncrna", "misc_RNA"),
                                            id=f.id, qualifiers={"locus_tag": f.id, "note": f.qualifiers["Note"]
                            , "description": f.qualifiers["description"]
                            , "gene": f.qualifiers["Note"]
                            , "alt name": f.qualifiers["description"]
                            , "product": f.qualifiers["description"]})
                elif f.type.lower() in ["contig", "exon", "cdsvi","CDS"]:
                    seqfeature = None
                else:
                    _log.warning("unknow feature " + f.type)
                    seqfeature = None
                if seqfeature:
                    record.features.append(seqfeature)
            return record

        with open(self.gb_file, "w") as h:
            for contig in tqdm(contig_iterator):
                if contig.features:
                    new_contig = process_contig(contig)
                    if new_contig.features:
                        bpio.write(new_contig, h, "gb")

    def execute(self, pwtools_path="pathway-tools", proxy="PROXY=proxy.fcen.uba.ar:8080"):
        """
         :param pwtools_path: complete path to pathway-tools binary. by default assumes that it is in the PATH
        """
        cmd = proxy + ' ' + pwtools_path + ' -no-cel-overview -no-web-cel-overview  -patho ' + self.workdir
        execute(cmd)

    def copy_qualifiers(self, mapping, f_in, f_out):
        for k, v in mapping.items():
            if k in f_in.qualifiers:
                if isinstance(v, dict):
                    for key, fn in v.items():
                        value = [y for y in [fn(x) for x in f_in.qualifiers[k]] if y]
                        if value:
                            f_out.qualifiers[key] = value

                elif isinstance(v, tuple):
                    dst = v[0]
                    value = [y for y in [v[1](x) for x in f_in.qualifiers[k]] if y]
                    if value:
                        f_out.qualifiers[dst] = value
                else:
                    dst = v
                    if f_in.qualifiers and (k in f_in.qualifiers):
                        value = [x for x in f_in.qualifiers[k]]
                        if value:
                            f_out.qualifiers[dst] = value

    def map_atributes(self, mapping, f_in, f_out):
        """
        "gene.mRna" : {
            "type" : "CDS",
            "qualifiers": {
                "gene_symbol" : "gene",
                "locus tag": "locus_tag",
                "x" : "db xref",
                "GO" : ("go_component", process_go_cc ),
        """
        for ftype, type_map in mapping.items():
            types = ftype.split(".")
            assert f_in.type == types[0]
            f_in1 = f_in
            if len(types) > 1 and hasattr(f_in, "sub_features") and f_in.sub_features:
                for subtype in types[1:]:
                    f_in1 = [x for x in f_in1.sub_features if x.type == subtype]
                    if f_in1:
                        f_in1 = f_in1[0]
                    else:
                        print (f_in)
                        continue
            if f_in1:
                for k, v in type_map.items():
                    if k == "qualifiers":
                        self.copy_qualifiers(v, f_in1, f_out)
                    else:
                        setattr(f_out, k, getattr(f_in, v) if v.startswith("$") else v)


if __name__ == "__main__":
    from SNDG import init_log
    from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB

    init_log()
    pw = PathwayTools("/data/organismos/GCF_001624625.1/annotation/")

    mdb = BioMongoDB("tdr")
    # contigmap = bpio.to_dict(
    #     bpio.parse("/data/organismos/cruzi/TriTrypDB-34_TcruziCLBrenerEsmeraldo-like_Genome.fasta", "fasta"))
    # pw.create_genebank(mdb.organism_iterator("cruzi", contigmap))
    # pw.create_organism_params(name="cruzi", organism="Trypanosoma cruzi strain CL Brener", domain="TAX-2759",
    #                           tax="353153")
    # pw.create_genetic_elements()
    # pw.execute("/opt/pathway-tools/pathway-tools")

    contigmap = {x.id.split(".")[0]: x for x in
                 bpio.parse("/data/organismos/GCF_001624625.1/annotation/GCF_001624625.1_ASM162462v1_genomic.gb", "gb")}
    pw.create_genebank(mdb.organism_iterator("GCF_001624625.1", contigmap))
    pw.create_organism_params(name="Bartonela", organism="Bartonella bacilliformis strain USM-LMMB 07", domain="TAX-2",
                              tax="774")
    pw.create_genetic_elements()
    pw.execute("/opt/pathway-tools/pathway-tools")
