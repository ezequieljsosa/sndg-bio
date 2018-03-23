'''
Created on Jul 20, 2017

@author: eze
'''

import json
import logging
import os
import shutil
import subprocess as sp

import Bio.SeqIO as bpio
from BCBio import GFF
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from SNDG import init_log, execute
from SNDG.BioMongo.Model.Sequence import Contig
from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
from bson.objectid import ObjectId

# from BIADeploy.BiaSql import BiaSql
# from SNDGInt.Submitter import ExternalAssembly


_log = logging.getLogger(__name__)


class JBrowse(object):
    '''
    classdocs
    '''

    def __init__(self, db, jbrowse_dir='/data/xomeq/jbrowse/',
                 path_PERL5LIB='/home/eze/perl5/lib/perl5'):
        '''
        Constructor
        '''
        self.jbrowse_dir = jbrowse_dir;
        self.db = db
        self.sequences = None
        self.path_PERL5LIB = path_PERL5LIB

    def load_sequences(self, path, seq_format="gb"):
        self.sequences = {x.id: x.seq for x in bpio.parse(path, seq_format)}

    def organism_dir(self, organism):
        return self.jbrowse_dir + "/data/" + organism + '/'

    def create_genome(self, organism, gff4jbrowse_fasta="/tmp/jbrowse_g.fasta", create_fasta=True,
                      extra_features={}):
        organism = organism.replace(' ', '_')

        if not os.path.exists(self.organism_dir(organism)):
            os.makedirs(self.organism_dir(organism))

        if create_fasta:
            with open(gff4jbrowse_fasta, "w") as h:
                p = {"name": 1}
                if not self.sequences:
                    p["seq"] = 1
                for i, c in enumerate(self.db.contig_collection.find({"organism": organism}, p)):
                    _log.debug("contig %i %s" % (i, c["name"]))
                    if self.sequences:
                        if c["name"] in self.sequences:
                            seq_id = c["name"]
                        else:
                            seq_id = [x for x in self.sequences if x.startswith(c["name"])][0]
                        bpio.write(SeqRecord(id=c["name"], seq=self.sequences[seq_id]), h, "fasta")
                    else:
                        if not c["seq"]:
                            raise Exception("Empty sequence")
                        bpio.write(SeqRecord(id=c["name"], seq=Seq(c["seq"])), h, "fasta")

        execute('PERL5LIB={perl5lib}  {base}bin/prepare-refseqs.pl --fasta "{fasta}" --out "{out_dir}" --key "{name}" ',
                perl5lib=self.path_PERL5LIB,
                fasta=gff4jbrowse_fasta, base=self.jbrowse_dir, out_dir=self.organism_dir(organism), name="Sequence")

        gff4jbrowse_file = "/tmp/jbrowse_g.gff"
        contigs = []
        proceced = {}
        for contig in Contig.objects(organism=organism).no_cache():

            seq_id = contig.name

            c = SeqRecord(id=seq_id, seq=Seq(""))
            if seq_id in extra_features:
                c.features = extra_features[seq_id]
            contigs.append(c)
            for f in contig.features:
                name = None
                if f._id:
                    ident = str(f._id)
                else:
                    ident = str(ObjectId())

                fl = FeatureLocation(start=f.location.start, end=f.location.end, strand=f.location.strand)
                ftype = f.type
                note = ""
                if f.type in ["gene", "CDS"]:
                    prot = self.db.proteins.find_one({"organism": organism, "alias": f.identifier},
                                                     {"gene": 1, "description": 1})
                    if not prot:
                        if ("tRNA" in f.identifier) or ("rRNA" in f.identifier) or ("ribosomal" in f.identifier):
                            name = f.identifier

                            ftype = "tRNA" if "tRNA" in f.identifier else "rRNA"
                        else:
                            _log.warn("gene %s was not found in the genome %s" % (f.identifier, organism))
                            continue
                    else:
                        name = " | ".join(list(set(prot["gene"][0:3])))
                        note = prot["description"]
                elif f.type in ["rRNA", "tRNA"]:
                    pass
                else:
                    print f

                if name in proceced:
                    continue
                proceced[name] = 1

                feature = SeqFeature(fl, type=ftype, id=ident,
                                     qualifiers={"ID": ident, "Name": name, "Note": note, "locus_tag": f.identifier})
                c.features.append(feature)

        with open(gff4jbrowse_file, "w") as h:
            GFF.write(contigs, h)

        os.chdir(self.jbrowse_dir)
        execute(
            'PERL5LIB=/home/eze/perl5/lib/perl5 ./bin/flatfile-to-json.pl --gff "{gff}" --out "{out_dir}" --key "{name}" --trackLabel "{name}" ',
            gff=gff4jbrowse_file, out_dir=self.organism_dir(organism), name="Genes")

        track_list_path = self.organism_dir(organism) + "/trackList.json"
        with open(track_list_path) as handle:
            txt = handle.read()

        str_on_click = '''
         "onClick"  : {
              "label": "go to product",
              "url": "function() { if(this.feature.get('type') == 'gene' ){  window.parent.location.href =  window.parent.location.href.split('genome/')[0] + '/protein/gene/' + this.feature.get('locus_tag') + '/'; }}"
          }
        '''
        txt = txt.replace('"key" : "Genes"', '"key" : "Genes",' + str_on_click)
        txt = txt.replace('"key":"Genes"', '"key" : "Genes",' + str_on_click)
        with open(track_list_path, "w") as handle:
            handle.write(txt)

        if not os.path.exists(self.organism_dir(organism) + "/names/"):
            os.makedirs(self.organism_dir(organism) + "/names/")
        with open(self.organism_dir(organism) + "/names/root.json", "w") as handle:
            handle.write('')

    #         sp.call('scp -r "' + self.organism_dir(organism) + '" 157.92.24.249:' + self.organism_dir(organism), shell=True)
    # --trackType CanvasFeatures

    def add_strain(self, organism, strain, vcf_path, bam_path):
        jb_vcf = self.organism_dir(organism) + strain + ".vcf"
        shutil.copy(vcf_path, jb_vcf)
        sp.call("bgzip -f " + jb_vcf, shell=True)
        sp.call("tabix -p vcf " + jb_vcf + ".gz", shell=True)

        jb_bam = self.organism_dir(organism) + strain + ".bam"
        shutil.copy(bam_path, jb_bam)
        sp.call("samtools index " + jb_bam, shell=True)

        track_list_path = self.organism_dir(organism) + "/trackList.json"
        shutil.copy(track_list_path, track_list_path + ".bk")
        with open(track_list_path) as handle:
            data = json.load(handle)
            data["tracks"].append(
                {
                    "label": strain + "_vcf",
                    "key": "SNPs " + strain,
                    "storeClass": "JBrowse/Store/SeqFeature/VCFTabix",
                    "urlTemplate": strain + ".vcf.gz",
                    "type": "JBrowse/View/Track/HTMLVariants"
                }

            )
            data["tracks"].append(
                {
                    "label": strain + "_bam",
                    "urlTemplate": strain + ".bam",
                    "type": "Alignments2"
                }
            )

        with open(track_list_path, "w") as handle:
            json.dump(data, handle,
                      indent=4, separators=(',', ': '))


if __name__ == "__main__":
    init_log()
    #     sqldb = BiaSql("biosql",user="root",passwd="mito")
    #     for i,x in enumerate(ExternalAssembly.select().where(ExternalAssembly.genome.is_null() )):
    #         print ( "procesando: " +  str(i))
    #         mdb = BioMongoDB("saureus")
    #
    #         if (not os.path.exists("/data/xomeq/jbrowse/data/" + x.assembly_accession + "/trackList.json") and
    #             mdb.seq_col_exists(x.assembly_accession)):
    #             jw = JBrowse(db=BioMongoDB("saureus").db)
    #             jw.sequences = {x.id:x.seq for x in sqldb.sequences(x.assembly_accession + "_ncbi")}
    #             jw.create_genome(x.assembly_accession)

    #     extra_features = defaultdict(lambda : [])
    #     for gff in open("/data/ger/ncbi_IP4.gff3.whole").read().split("##gff-version 3")[1:]:
    #         for r in GFF.parse(StringIO.StringIO(gff)):
    #             for f in r.features:
    #                 if f.type == "TATA":
    #                     extra_features[r.id].append(f)
    #                 elif f.type == "TSS":
    #                     extra_features[r.id].append(f)
    mdb = BioMongoDB("saureus")
    jw = JBrowse(db=mdb.db)
    for genome in ["HIV-1","EcoliMG1655","SacCereS288C"]:
        if not os.path.exists("/data/xomeq/jbrowse/data/"+ genome + "/xomeq" ):
            try:
                jw.create_genome(genome)
            except Exception as ex:
                _log.warn(ex)
    # jw.sequences = {x.id:x.seq for x in bpio.parse("/data/projects/Staphylococcus/annotation/ncbi/GCF_000009645.1_ASM964v1_genomic.gb","gb")}
    # jw.create_genome("SAureusN315") #gff4jbrowse_fasta="/data/organismos/ILEX_PARA/xomeq/jbrowse_g.fasta", create_fasta=False, extra_features=dict(extra_features)


    # jw.load_sequences("/data/organismos/cruzi/TriTrypDB-34_TcruziCLBrenerEsmeraldo-like_Genome.fasta",seq_format="fasta")
    # jw.create_genome("cruzi")
    #
    # jw.load_sequences("/data/organismos/Pext14-3B/annotation//GCF_000242115.1_Pext14-3B_1.0_genomic.gbff")
    # jw.create_genome("Pext14-3B")
    #
    # jw.load_sequences("/data/organismos/GCF_001624625.1/annotation/GCF_001624625.1_ASM162462v1_genomic.gb")
    # jw.create_genome("GCF_001624625.1")

#     for s in [ "15-6324_S3_L001","2003_S4_L001"]:
#         vcf = "/data/projects/PiuriTB/analysis/variant_call_h37/" + s + "/variants.ann.vcf"
#         bam = "/data/projects/PiuriTB/analysis/reads_h37rv_aln/" + s + "/final_bwa.bam"
#         jw.add_strain("H37Rv",s, vcf , bam)
