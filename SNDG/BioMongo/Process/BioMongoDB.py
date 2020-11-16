'''
Created on Apr 27, 2017

@author: eze
'''
import glob
import logging
import os
from itertools import combinations
from tqdm import tqdm

import numpy as np
import pymongo
from bson.objectid import ObjectId
from mongoengine import connect, register_connection
from tqdm import tqdm

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import Bio.SeqIO as bpio
from SNDG import init_log
from SNDG.BioMongo.Model.Protein import Protein
from SNDG.BioMongo.Model.Feature import Feature,Location
from SNDG.BioMongo.Model.SeqColDruggabilityParam import SeqColDruggabilityParam
from SNDG.BioMongo.Model.SeqCollection import SeqCollection, DataUpload
from SNDG.BioMongo.Model.Sequence import BioProperty, Contig
from SNDG.BioMongo.Process.BioCyc2Mongo import BioCyc
from SNDG.BioMongo.Process.COG2Mongo import COG2Mongo
from SNDG.BioMongo.Process.EC2Mongo import EC2Mongo
from SNDG.BioMongo.Process.GO2Mongo import GO2Mongo
from SNDG.BioMongo.Process.PFAM2Mongo import PFAM2Mongo

_log = logging.getLogger(__name__)


class FilesystemResolver():
    def __init__(self, basefs='/data/'):
        self.basefs = basefs

    def pathways_dir(self, name):
        pwdir = self.basefs + "/organismos/" + name + "/annotation/pathways/"
        assert os.path.exists(pwdir), pwdir
        return pwdir


class BioMongoDB(object):
    demo_id = ObjectId("563b9440b1b50423d1fd1fee")
    demo = "demo"
    GENE_FIELD_IMPORT = "id"

    def __init__(self, dbname, port=27017, host='127.0.0.1', basefs='/data/'):

        self.db = pymongo.MongoClient(host=host, port=port)[dbname]
        self.struct_db = pymongo.MongoClient(host=host, port=port)["pdb"]
        self.fs_resolver = FilesystemResolver(basefs)
        connect(dbname, host=host, port=port)
        register_connection("pdb", name="pdb", host=host, port=port)

        self.paths = {
            "gene": []
        }

    def seq_col_exists(self, seq_col_name):
        return SeqCollection.objects(name=seq_col_name).count()

    def delete_seq_collection(self, name):
        genome = self.db.sequence_collection.find_one({"name": name})
        if genome:
            genome_id = genome["_id"]
            print (self.db.sequence_collection.remove({"name": name}))
            print (self.db.contig_collection.remove({"seq_collection_id": genome_id}))
            print (self.db.structures.remove({"organism": name}))
            print (self.db.proteins.remove({"organism": name}))
            print (self.db.col_ont_idx.remove({"seq_collection_name": name}))
            print (self.db.users.update({}, {"$pull": {"links": {"description": name}}}, multi=True))
            print (self.db.users.update({}, {"$pull": {"seq_collections": genome_id}}, multi=True))

    def init_db(self, ec=True, pfam=True, tigrfam=True, cog=True, so=True, pathways=True, go=True):
        if ec:
            self.db.ontologies.remove({"ontology": "ec"})
            _log.debug("Loading EC...")
            ec = EC2Mongo(self.db)
            ec.load_enzclass("/data/databases/ec/enzclass.txt")
            ec.load_enzdata("/data/databases/ec/enzyme.dat")
            ec.load_children()
            _log.debug("EC loaded")

        if pfam:
            self.db.ontologies.remove({"ontology": "pfam"})
            _log.debug("Loading PFAM;...")
            pfam = PFAM2Mongo("/data/databases/xfam/Pfam-A.hmm", "pfam")
            pfam.load()
            _log.debug("PFAM loaded")

        if tigrfam:
            self.db.ontologies.remove({"ontology": "tigrfam"})
            _log.debug("Loading TigrFAM;...")
            pfam = PFAM2Mongo("/data/databases/xfam/TIGRFAMs_15.0_HMM.LIB", "tigrfam")
            pfam.load()
            _log.debug("TigrFAM loaded")

        if cog:
            self.db.ontologies.remove({"ontology": "cog"})
            cog = COG2Mongo()
            cog.load_fun("/data/databases/cog/fun.txt")
            cog.load_whog("/data/databases/cog/whog")

        if so:
            self.db.ontologies.remove({"ontology": "so"})
            _log.debug("Loading SO")
            go = GO2Mongo("/data/databases/so/so-xp.obo", self.db, ontology_name="so")
            go.load(False)
            _log.debug("SO loaded")

        if pathways:
            self.db.ontologies.remove({"ontology": "biocyc_pw"})
            self.db.ontologies.remove({"ontology": "biocyc_reac"})
            biocyc = BioCyc(self.db)
            _log.debug("Loading Biocyc...")

            for x in glob.glob("/data/databases/biocyc/*cyc/"):
                _log.debug(x)
                biocyc.load_pathways(x + "/pathways.col", database=x.split(os.sep)[-2])
                biocyc.load_reactions(x + "/reactions.dat", database=x.split(os.sep)[-2])
                biocyc.load_compounds(x + "/compounds.dat", database=x.split(os.sep)[-2])

            _log.debug("Biocyc loaded")
        if go:
            self.db.ontologies.remove({"ontology": "go"})
            _log.debug("Loading GO")
            go = GO2Mongo("/data/databases/go/go.obo", self.db, slim_file="/data/databases/go/goslim_generic.obo")
            go.load(False)
            _log.debug("GO loaded")

    def copy_genome(self, name, newname=None, dst_db=None):
        """
        :param name: of the collection in the source db
        :param newname: of the collection
        :param dst_db: destination database, by default is the same as the source db

        """

        sc = self.db.sequence_collection.find_one({"name": name})
        genome_id = sc["_id"]
        new_id = genome_id if dst_db else ObjectId()

        assert bool(newname) | bool(dst_db), "newname and db cant be none at the same time"

        if not dst_db:
            dst_db = self.db


        if not newname:
            newname = name

        sc["name"] = newname
        dst_db.sequence_collection.insert(sc)

        print("Copiando contigs:")
        total = self.db.contig_collection.count({"seq_collection_id": genome_id})
        for contig in tqdm(self.db.contig_collection.find({"seq_collection_id": genome_id}), total=total):
            contig["seq_collection_id"] = new_id
            contig["organism"] = newname
            dst_db.contig_collection.save(contig)

        print("Copiando proteinas:")
        total = self.db.proteins.count({"organism": name})
        for protein in tqdm(self.db.proteins.find({"organism": name}), total=total):
            protein["seq_collection_id"] = new_id
            protein["seq_collection_name"] = newname
            dst_db.proteins.save(protein)

        print("Copiando indices:")
        total = self.db.col_ont_idx.count({"seq_collection_name": name})
        for idx in tqdm(self.db.col_ont_idx.find({"seq_collection_name": name}), total=total):
            idx["seq_collection_id"] = new_id
            idx["seq_collection_name"] = newname
            dst_db.col_ont_idx.save(idx)

        print("Copiando estructuras:")
        total = self.struct_db.structures.count({"organism": name})
        for struct in tqdm(self.struct_db.structures.find({"organism": name}), total=total):
            dst_db.client.pdb.structure.save(struct)



    def delete_feature_type(self, organism, feature_type):
        self.db.proteins.update({"organism": organism, "features.type": feature_type},
                                {"$pull": {"features": {"type": feature_type}}},
                                multi=True)

    #         for p in Protein.objects(organism=organism).no_cache().timeout(False):
    #             p.features = [x for x in p.features if x.type != feature_type]
    #             p.save()

    def properties_from_feature(self, organism, feature_type, value_fn, property_name=None, url_fn=None):
        proteins = Protein.objects(__raw__={"organism": organism, "features.type": feature_type}).no_cache().timeout(
            False)
        if property_name == "human_offtarget":
            res = self.db.proteins.update({"organism": organism, "search.human_offtarget": {"$exists": True}},
                                          {"$unset": {"search.human_offtarget": ""}}, multi=True)
            _log.info(res)
        for p in proteins:
            f = [x for x in p.features if x.type == feature_type]
            if f:
                f = f[0]
                bp = BioProperty(_type=feature_type,
                                 property=property_name if property_name else feature_type,
                                 value=value_fn(f),

                                 )
                if url_fn:
                    bp.url = url_fn(f)
                props = [x for x in p.properties if
                         not (x._type == feature_type and ((not property_name) or (x.property == property_name)))]
                p.properties = props
                p.properties.append(bp)
                if property_name == "human_offtarget":
                    p.search.human_offtarget = 1 - f.aln.identity
                p.save()
        if property_name == "human_offtarget":
            res = self.db.proteins.update({"organism": organism, "search.human_offtarget": {"$exists": False}},
                                          {"$set": {"search.human_offtarget": 1}}, multi=True)
            _log.info(res)

    def clean_structures(self, organism):
        proteins = list(Protein.objects(organism=organism).no_cache().timeout(False))
        for p in proteins:
            experimentals = p.cristals()
            models = p.models()
            if experimentals and models:
                for model in models:
                    for exp in experimentals:
                        hit = [f for f in p.features if f.identifier.startswith(exp.name)][0]
                        model_query = model.templates[0].aln_query
                        model_range = set(range(model_query.start, model_query.end))
                        exp_range = set(range(hit.location.start, hit.location.end))
                        if ((len(model_range & exp_range) * 1.0 / len(model_range)) > 0.8):
                            model.delete()

            if len(models) > 1:
                for i, j in combinations(range(len(models), 2)):
                    if ((models[i].templates[0].aln_query.start == models[j].templates[0].aln_query.start)
                            and (models[i].templates[0].aln_query.end == models[j].templates[0].aln_query.end)):
                        models[i].delete()

    def props_from_dbxref(self, name):

        i = 0
        for p in Protein.objects(organism=name).no_cache().timeout(False):
            if p.dbxrefs:
                i += 1
                prop = BioProperty(_type="dbxref", property="links", value=p.dbxrefs)
                p.alias += [x.split(":")[1] for x in p.dbxrefs if x.lower().startswith("uniprot")]
                p.alias = list(set(p.alias))
                p.properties.append(prop)
                p.save()


    def load_metadata(self, organism_name, datafile, uploader=demo):
        import pandas as pd
        from tqdm import tqdm

        seqCollection = list(SeqCollection.objects(name=organism_name))
        seqCollection = seqCollection[0]
        errors = []

        upload = DataUpload(uploader=uploader, errors=errors);

        df = pd.read_table(datafile, comment="#", index_col=False)

        headerProperties = [c for c in df.columns if c != BioMongoDB.GENE_FIELD_IMPORT]
        prots = Protein.objects(organism=organism_name)
        for hp in headerProperties:
            prots.update(
                __raw__={"$pull": {"properties": {"property": hp, "_type": uploader}}, "$unset": {"search." + hp: ""}})

        upload.properties = headerProperties

        numericFields = []

        for k, v in dict(df.dtypes).items():
            if v not in [np.float64, np.int64]:
                df[k] = df[k].astype('category')
            else:
                numericFields.append(k)

        assert BioMongoDB.GENE_FIELD_IMPORT in df.columns

        for linenum, fields in tqdm(df.iterrows()):

            gene = fields[BioMongoDB.GENE_FIELD_IMPORT]

            if not gene:
                text = str(linenum) + " gene field is empty"
                errors.append(text);
                continue

            count = Protein.objects(organism=organism_name, alias=gene).count();

            if not count:
                text = str(linenum) + " " + gene + " does not exists in " + organism_name
                print(text)
                errors.append(text);
                continue;

            prots = Protein.objects(organism=organism_name, alias=gene)

            for propertyName in headerProperties:
                prop = {"_type": uploader, "value": fields[propertyName]}
                prop["property"] = propertyName
                prots.update(__raw__={
                    "$push": {"properties": prop},
                    "$set": {"search." + propertyName: fields[propertyName]}}
                )

        for p in headerProperties:
            dpType = "number" if p in numericFields else "value"

            options = [] if p in numericFields else list(set(df[p]))
            currentDp = seqCollection.druggabilityParam(p, uploader)

            if currentDp:
                currentDp = currentDp[0]
                currentDp.options = options
                currentDp.type = dpType
            else:
                dp = SeqColDruggabilityParam(type=dpType, name=p, options=options,
                                             uploader=uploader, target="protein")
                seqCollection.druggabilityParams.append(dp)

        seqCollection.uploads.append(upload);
        seqCollection.save()

    @staticmethod
    def protein_fasta(outfile_path, organism):
        with open(outfile_path, "w") as h:
            for p in Protein.objects(organism=organism).no_cache():
                r = SeqRecord(id=p.gene[0], description="", seq=Seq(p.seq))
                bpio.write(r, h, "fasta")

    def organism_iterator(self, organism, seq_map=None):
        for dbcontig in Contig.objects(organism=organism).no_cache():
            if seq_map:
                seq = str(seq_map[dbcontig.name].seq)
            else:
                seq = dbcontig.seq
            contig = SeqRecord(id=dbcontig.name, seq=Seq(seq))
            for dbfeature in dbcontig.features:
                qualifiers = {"locus_tag": [dbfeature.locus_tag]}
                p = list(Protein.objects(organism=organism, gene=dbfeature.identifier))
                if p:
                    p = p[0]
                    qualifiers["description"] = [p.description]
                    qualifiers["gene_symbol"] = p.gene
                    qualifiers["Note"] = [p.description]

                    ecs = [x.upper() for x in p.ontologies if x.startswith("ec:")]
                    gos = [x.upper() for x in p.ontologies if x.startswith("go:")]
                    if ecs:
                        qualifiers["EC"] = ecs
                    if gos:
                        qualifiers["GO"] = gos
                    feature = SeqFeature(id=dbfeature.identifier, type=dbfeature.type, qualifiers=qualifiers,
                                         location=FeatureLocation(start=dbfeature.location.start,
                                                                  end=dbfeature.location.end,
                                                                  strand=dbfeature.location.strand))
                    contig.features.append(feature)
            yield contig

    def load_from_emapper(self, organism, emapperv2_file):
        from SNDG.Annotation.EMapper import EMapper
        em = EMapper()
        em.read_file(emapperv2_file)
        for locus_tag, record in em.data.items():
            prot = Protein.objects(organism=organism, gene=locus_tag).get()
            for ec in record["EC"].split(","):
                prot.ontologies.append("ec:" + ec)
            for go in record["GOs"].split(","):
                prot.ontologies.append(go.lower())
            prot.save()

    def load_from_interpro(self, organism, interprot_gff):
        for l in tqdm(open(interprot_gff)):
            if l.startswith(">"):
                break
            if l.startswith("##"):
                continue
            l = l.replace("EC=", "EC ")
            locus_tag, source, feature, start, end, score, strand, frame = l.split("\t")[:8]
            attributes = " ".join(l.split("\t")[8:])

            if feature == "polypeptide":
                continue

            start, end = int(start), int(end)

            if "signature_desc=" in attributes:
                repl = attributes.split("signature_desc=")[1].split(";Name=")[0]
                attributes = attributes.replace(repl,repl.replace("=","%3D").replace(";","%3B"))

            attributes = {x.split("=")[0]: x.split("=")[1] for x in attributes.split(";")}
            # [seq,source,feature,start,end,score,strand,frame,attributes ])
            feature = Feature(_id=ObjectId(), location=Location(start=start, end=end),
                              identifier=attributes["Name"], type=source)
            prot = Protein.objects(organism=organism, gene=locus_tag).get()



            if "signature_desc" in attributes:
                feature.qualifiers = {"description":attributes["signature_desc"]}
            if "Ontology_term" in attributes:
                for ont in attributes["Ontology_term"].split(","):
                    ont = ont.replace('"',"").strip()
                    prot.ontologies.append(ont.lower())
            if "Dbxref" in attributes:
                for ont in attributes["Dbxref"].split(","):
                    ont = ont.replace('"',"").strip()
                    prot.ontologies.append(ont.lower())

            prot.features.append(feature)
            prot.save()

import re

regx = re.compile("^ec:", re.IGNORECASE)

if __name__ == '__main__':
    init_log()
    # pdb = pymongo.MongoClient().pdb

    mdb = BioMongoDB("tdr")

    mdb.load_metadata("Kp13", "/home/eze/workspace/kp13/metadata.tbl")

    # tax_db.initialize(MySQLDatabase('bioseq', user='root', passwd="mito"))
    #     bacs = []
    #     ts = []
    #     for x in open("/data/databases/biocyc/metacyc/pathways.dat"):
    #         if x.strip().startswith("TAXONOMIC-RANGE"):
    #             try:
    #                 t = x.split("TAX-")[1].strip()
    #                 ts.append(t)
    #             except IndexError:
    #                 pass
    #     ts = set(ts)
    #     for t in ts:
    #         try:
    #             tax = Tax.getTax(t)
    #             if tax:
    #                 try:
    #                     root = Tax.parents(tax)[-1]
    #                     if root.ncbi_taxon_id == 2:
    #                         bacs.append(t)
    #                 except IndexError:
    #                     pass
    #         except :
    #             pass
    #     print bacs

    # from BIADeploy.BioDocFactory import BioDocFactory
    # fasta_path = "/data/ger/ncbi_LBH.fna2"
    # gff_path = "/data/ger/ncbi_LBH.gff32"
    # mdb.delete_seq_collection("SAureusN315")
    #     mdb.load_genome_fasta_and_gff("L_HELV", fasta_path, gff_path)

    #     mdb.index_seq_collection("ILEX_PARA",pathways=False,go=False,keywords=True,ec=True,organism_idx=True)

    gff_path = "/data/ger/ncbi_IP4.gff3.whole"

    #     bcyc = BioCyc("xomeq")
    #     bcyc.pre_build_index(SeqCollection.objects(name="ILEX_PARA").get())

    #     for i,sq in enumerate(mdb.db.sequence_collection.find({},{"name":1,"organism":1,"statistics":1,"ncbi_assembly":1})):
    #         stats = len(sq["statistics"]) if "statistics" in sq else 0
    #         if not stats:
    #             prots = mdb.db.proteins.count({"organism":sq["name"],"keywords.0":{"$exists":1}})
    #             print  sq["ncbi_assembly"] if "ncbi_assembly" in sq else sq["name"] + " " + str(prots)

    #     pepe = 0
    #     for s in pdb.structures.find({"organism":{"$exists":1},"tax":{"$exists":0}},{"organism":1},no_cursor_timeout=True):
    #         org = s["organism"].replace(":","").replace(";","").replace("SYNTHETIC CONSTRUCT","").strip()
    #         if org:
    #             for tn in TaxName.select( ).where((fn.Lower(TaxName.name) ** ("%" + org.lower() + "%")) & (TaxName.name_class == TaxName.scientific_name)):
    #                 pepe = pepe + 1
    #                 try:
    #                     pdb.structures.update({"_id":s["_id"]},{"$set":{"tax":int(tn.taxon.ncbi_taxon_id)}})
    #                 except:
    #                     pass
    #                 break
    #     print pepe

    #     mdb.load_metadata("H37Rv","/media/eze/Data/data/projects/patho/H37Rv/annotation/vfdb.tbl")

    #     collection = SeqCollection.objects(name="Pext14-3B").get()
    #     indexer = StructuromeIndexer(collection)
    #     indexer.build_index()

    #     mdb.delete_seq_collection("LHelv")
    #     mdb.load_from_db("LHelv", "LHelv",NCBI.f_mRNA)
    #     mdb.load_hmm("Pext14-3B","/data/organismos/Pext14-3B/annotation/pfam/domains.hmm")
    #     mdb.load_pdb_domains("Pext14-3B", "/data/organismos/Pext14-3B/annotation/pdb/dns_from_pdb.xml")
    #     mdb.load_pathways("Pext14-3B","/data/organismos/Pext14-3B/annotation/pathways/pathways.sbml")
    #     mdb.index_seq_collection("Pext14-3B")
    #     mdb.build_statistics("Pext14-3B")

    # mdb.index_seq_collection(genome, ec, go, keywords, organism_idx, pathways)

    #     sa = StructureAnotator()
    #
    #     for i,total,model in sa.iterator(pymongo.MongoClient().pdb, "TGONDII",{"name":"41.m00001.1"}):
    #
    #         try:
    #             _log.debug( model.name + ":" +  str(i) + "/" + str(total))
    #             template = model.templates[0]
    #             protein = Protein.objects(organism="TGONDII",alias=template.aln_query.name).get()
    #             sa.annotate_model(model,protein.domains())
    #             model.save()
    #
    #         except Exception:
    #             import traceback
    #             traceback.print_exc()

    #     tbd = TBDream()
    #     tbd.load()
    #     tbd.load_in_sndg()

    # mdb.index_seq_collection("H37Rv", ec=True, go=True, keywords=False, organism_idx=False, pathways=False)
    # mdb.build_statistics("H37Rv")
    #     mdb.delete_seq_collection("ILEX_PARA")
    #     mdb.load_genome_fasta_and_gff("ILEX_PARA", "/data/ncbi_IP3_2.fna", "/data/ncbi_IP3_2.gff3")

    #     for strain in ["Ra"]: #,
    # #                 "0037","0058","0142","0271","0298","0450","0484","0564","1096","1300",
    # #                  "1445","1493","1527","1584","1707","1710","1764",
    # #                  "1796","1875","3296","3867INF","3867NE","3867NI"]:
    # #         mdb.annotate_variants("H37Rv", strain, "tbdream",TBDream.parse_change)
    #         ca = CoverageAnalysis()
    # #         ca.run_sdepth("/data/projects/PiuriTB/analysis/reads_h37rv_aln/" + strain + "/final_bwa.bam")
    # #         ca.min_depth = 5
    # #         ca.gb_cov("/data/organismos/H37Rv/annotation/ncbi/h37rv.3.gb",fname=lambda f:f.qualifiers["locus_tag"][0] if "locus_tag" in f.qualifiers else f.id)
    #
    #         prop = strain + "_low_cov_variant"
    #         collection = SeqCollection.objects(name="H37Rv").get()
    #         Protein.objects(organism="H37Rv").update(__raw__={"$set":{"search." + prop :False}})
    #         if not collection.has_druggability_param(prop):
    #                 dp = SeqColDruggabilityParam(name=prop, description="There is a reported variant with low coverage in this strain",
    #                                              target="variant-strain",  type=SeqColDruggabilityParamTypes.value, uploader="demo")
    #                 dp.options = ["true", "false"]
    #                 dp.defaultValue = "true"
    #                 dp.defaultOperation = "equal"
    #                 dp.defaultGroupOperation = "avg"
    #                 collection.druggabilityParams.append(dp)
    #         collection.save()

    #         for p in Protein.objects(__raw__={"organism":"H37Rv","search.resistance":True}).no_cache():
    #             for f in p.features:
    #                 if f.type == "tbdream":
    #                     rows = ca.depth_df[ca.depth_df.pos == f.location.start]
    #
    #                     if  rows.empty or rows.iloc[0].depth < 20:
    #                         Protein.objects(organism="H37Rv",gene=p.gene[0]).update(__raw__={"$set":{"search." + prop :True}})

    #         prop = strain + "_low_cov"
    #         collection = SeqCollection.objects(name="H37Rv").get()
    #         Protein.objects(organism="H37Rv").update(__raw__={"$set":{"search." + prop :False}})
    #         if not collection.has_druggability_param(prop):
    #                 dp = SeqColDruggabilityParam(name=prop, description="The protein has low coverage with " + strain + " reads mapping",
    #                                              target="variant-strain",  type=SeqColDruggabilityParamTypes.value, uploader="demo")
    #                 dp.options = ["true", "false"]
    #                 dp.defaultValue = "true"
    #                 dp.defaultOperation = "equal"
    #                 dp.defaultGroupOperation = "avg"
    #                 collection.druggabilityParams.append(dp)
    #         collection.save()
    #         ca.feature_cov_df = pd.read_csv("/data/projects/PiuriTB/analysis/reads_h37rv_aln/" + strain + "/coverage.csv",index_col=False)
    #         for f in ca.not_cov_features(0.75):
    #             Protein.objects(organism="H37Rv",gene=f.id).update(__raw__={"$set":{"search." + prop :True}})

    #
    #     mdb.annotate_variants_with_prots("SaureusN315",["Aanensen2016"],Saureus.drugs,True)
    #     mdb.indexVariants("SaureusN315")

    # mdb.annotate_variants_with_prots("H37Rv",["tbdream"],TBDream.drugs,True)
    # mdb.indexVariants("H37Rv")

    #     for x in ["TB-16-1109","TB-15-2748","TB-15-6324","TB-16-138"]:
    #         mdb.delete_seq_collection(x)


