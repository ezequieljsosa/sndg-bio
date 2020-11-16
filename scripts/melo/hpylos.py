import logging
from SNDG import Struct, init_log, mkdir

init_log()
from glob import glob
import pymongo
from SNDG.WebServices.NCBI import NCBI
from SNDG.WebServices.Offtargeting import Offtargeting
import Bio.SeqIO as bpio
from tqdm import tqdm
import json
import sys
logging.getLogger("peewee").setLevel(logging.WARN)
from peewee import MySQLDatabase
from SNDG.BioMongo.Process.Taxon import tax_db
from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
from SNDG.Sequence.ProteinAnnotator import ProteinAnnotator, Mapping
from SNDG.BioMongo.Process.Importer import from_ref_seq, update_proteins, import_prop_blast
from SNDG.BioMongo.Process.BioDocFactory import BioDocFactory
from SNDG.BioMongo.Model.Protein import Protein,ChEMBL
from SNDG.Network.KEGG import Kegg
from SNDG.BioMongo.Process.Importer import _common_annotations, _protein_iter, import_kegg_annotation, \
    index_seq_collection, build_statistics, load_pathways
from BCBio import GFF
from SNDG.BioMongo.Process.Taxon import Tax
from SNDG.BioMongo.Model.Structure import ModeledStructure, Molecule, ResidueAln, SimpleAlignment, StructureQuality, \
    ExperimentalStructure, Chain,SeqCollection
from SNDG.BioMongo.Model.Alignment import AlnLine
import os
from SNDG.BioMongo.Process.StructureAnotator import StructureAnotator
import Bio.SearchIO as bpsio
from Bio.SeqUtils import seq1,seq3
tax_db.initialize(MySQLDatabase('bioseqdb', user='root', passwd="mito"))
mdb = BioMongoDB("tdr", port=27017)
mysqldb = ProteinAnnotator.connect_to_db(database="unipmap", user="root", password="mito")

orgs = [("Mpylori26695", "Helicobacter pylori 26695 (e-proteobacteria)",
         "/data/organismos/Mpylori26695/GCF_000008525.1_ASM852v1_genomic.gbff", 85962),
        ("MpyloriIndia", "Helicobacter pylori India7 (e-proteobacteria)",
         "/data/organismos/MpyloriIndia/GCF_000185185.1_ASM18518v1_genomic.gbff", 907238),

        ]

for name, org, ann_path, tax in orgs:
    organism = name
    mkdir("/data/organismos/" + name + "/annotation/offtarget")
    mkdir("/data/organismos/" + name + "/annotation/pwtools")
    mkdir("/data/organismos/" + name + "/annotation/pathways")
    mkdir("/data/organismos/" + name + "/estructura/raw")
    mkdir("/data/organismos/" + name + "/estructura/sndg/modelos")
    mkdir("/data/organismos/" + name + "/estructura/sndg/pockets")

    from_ref_seq(name, ann_path, tax=tax, cpus=3)
    mdb.protein_fasta("/data/organismos/" + name + "/annotation/proteins.faa", name)
    update_proteins("/tmp/" + name + "/", "/data/organismos/" + name + "/annotation/proteins.faa", name, 1003200,
                    db_init=mysqldb)

    Offtargeting.offtargets("/data/organismos/" + name + "/annotation/proteins.faa",
                            "/data/organismos/" + name + "/annotation/offtarget/",
                            offtarget_dbs=["/data/databases/deg/degaa-p.dat",
                                           "/data/databases/human/gencode.v17.pc_translations.fa",
                                           "/data/databases/human/gut_microbiota.fasta"]
                            )
    import_prop_blast(mdb.db, name, "hit_in_deg",
                      "/data/organismos/" + name + "/annotation/offtarget/degaa-p.tbl",
                      "table", "Hit in DEG database",
                      value_fn=lambda x: x.identity > 70,
                      default_value="true",
                      no_hit_value="false", choices=["true", "false"], type="value", defaultOperation="equal")

    import_prop_blast(mdb.db, name, "human_offtarget",
                      "/data/organismos/" + name + "/annotation/offtarget/gencode.tbl",
                      "table", "Human offtarget score (1 - best hit identity)",
                      value_fn=lambda x: 1 - (x.identity * 1.0 / 100),
                      default_value=0.4,
                      no_hit_value=1)
    import_prop_blast(mdb.db, name, "gut_microbiota_offtarget",
                      "/data/organismos/" + name + "/annotation/offtarget/gut_microbiota.tbl",
                      "table", "Gut microbiota offtarget score (1 - best hit identity)",
                      value_fn=lambda x: 1 - (x.identity * 1.0 / 100),
                      default_value=0.4,
                      no_hit_value=1)

from Bio.SeqUtils import seq1
from Bio.PDB.PDBParser import PDBParser
import logging
import re

_log = logging.getLogger(__name__)
with open('/data/databases/pdb/manually_curated/compound_type.csv') as handle:
    compound_type = {x.replace('"', "").strip().split(",")[0]: x.replace('"', "").strip().split(",")[1] if
    x.replace('"', "").strip().split(",")[1] else "?" for x in handle.readlines()}


def get_compound_type(residue):
    return compound_type[residue.get_resname().strip()] if residue.get_resname().strip() in compound_type else "?"


from SNDG.Structure.FPocket import FPocket


def hunt_pockets(pdb):
    pocket_data = pdb + ".pocket.json"
    if not os.path.exists(pocket_data):
        fpo = FPocket(pdb)
        result = fpo.hunt_pockets()
        result.save(pocket_data)
        result.delete_dir()


def build_assessments(pdb_dir, pdb_name):
    cmd = """docker run --rm  -v %s:/out structurome  modellerkey python -c 'import json;from SNDG.Structure.QMean import QMean;assessment = QMean.assesment("/out/%s");
with open( "/out/%s.json", "w") as h:
    json.dump(assessment, h)'""" % (pdb_dir, pdb_name, pdb_name)
    sp.check_output(cmd, shell=True)


parser = PDBParser(PERMISSIVE=1, QUIET=1)
import subprocess as sp

for name, org, ann_path, tax in orgs:
    organism = name
    basepath = "/data/organismos/" + name + "/estructura/raw/"
    models_dir = "/data/organismos/" + name + "/estructura/sndg/modelos/"
    for x in glob(basepath + "/*/*/*/*.pdb"):
        model_file = models_dir + x.split("/")[-2] + ".pdb"
        if not os.path.exists(model_file):
            sp.call(["ln", x, model_file])

    model_files = glob(models_dir + "*.pdb")
    seq_col_id = SeqCollection.objects.get(name=organism).id
    with tqdm(model_files) as pbar:
        for model_file in pbar:
            pbar.set_description("processing %s" % model_file)

            model_filename = model_file.split("/")[-1]
            model_name = model_filename.split(".pdb")[0]

            hunt_pockets(model_file)
            build_assessments(os.path.dirname(model_file), os.path.basename(model_filename))

            seq_name = "_".join(model_name.split("_")[:-4])
            org_model_name = model_name
            template_name = "_".join(model_name.split("_")[-4:])

            if ModeledStructure.objects(organism=organism, name=model_name).count():
                continue

            prot = list(Protein.objects(organism=organism, gene=seq_name))
            if len(prot) == 0:
                _log.warn("Not found: " + seq_name)
                continue

            aln = [hit[0] for hit in list(bpsio.read(basepath + "/" + seq_name + "/profile_search.xml", "blast-xml")) if
                   hit.id == template_name][0]

            with open(model_file + ".json") as h:
                assessments = json.load(h)
            pockets = []

            try:
                structure = parser.get_structure(model_name, model_file)
                with open(model_file) as h:
                    pdblines = [l.strip() for l in h.readlines() if l.startswith("REMARK")]
                strdoc = ModeledStructure(name=model_name, seq_collection_id=seq_col_id, pipeline="sndg",
                                          organism=organism)

                squality = StructureQuality(name="ga341", value=float(
                    [l.split("GA341 score: ")[1] for l in pdblines if "GA341 score: " in l][0]))
                strdoc.qualities.append(squality)
                squality = StructureQuality(name="dope", value=float(
                    [l.split("DOPE score: ")[1] for l in pdblines if "DOPE score: " in l][0]))
                strdoc.qualities.append(squality)
                squality = StructureQuality(name="qmean", value=assessments["QMEAN4_norm"])
                strdoc.qualities.append(squality)
                squality = StructureQuality(name="zqmean", value=assessments["QMEAN4_zscore"])
                strdoc.qualities.append(squality)

                model = structure[0]
                for chain in model:
                    chaindoc = Chain(name=chain.id)
                    strdoc.chains.append(chaindoc)

                    aln_h = str(aln.aln[1].seq)
                    aln_q = str(aln.aln[0].seq)
                    prot_start = aln.query_start
                    prot_end = aln.query_end
                    pdb_txt = [l for l in pdblines if "TEMPLATE: " in l][0].split()
                    pdb_start = int(pdb_txt[4].split(":")[0])
                    pdb_end = int(pdb_txt[6].split(":")[0])

                    hit_start = aln.hit_start
                    hit_end = aln.hit_end

                    aln_query = AlnLine(name=seq_name, seq=aln_q, start=prot_start, end=prot_end)
                    aln_hit = AlnLine(name=template_name, seq=aln_h, start=hit_start, end=hit_end)
                    template = ResidueAln(aln_query=aln_query, aln_hit=aln_hit,
                                          query_res_start=1,
                                          query_res_end=len(chain),
                                          hit_res_start=pdb_start,
                                          hit_res_end=pdb_end)
                    strdoc.templates.append(template)

                    for residue in chain:
                        molecule = Molecule(resid=residue.id[1],
                                            chain=chain.id,
                                            compound=seq1(residue.get_resname()),
                                            compound_type=get_compound_type(residue))
                        if molecule.compound_type == 'RESIDUE':
                            chaindoc.residues.append(molecule)
                        else:
                            if not [x for x in strdoc.ligands if
                                    (
                                            x.compound_type == molecule.compound_type) and molecule.compound_type == 'SOLVENT']:
                                strdoc.ligands.append(molecule)

                pockets_json = model_file + ".pocket.json"
                if os.path.exists(pockets_json):
                    rss = StructureAnotator.pocket_residue_set(pockets_json, model.get_atoms())
                    strdoc.pockets = rss
                strdoc.save()
            except ImportError as ex:
                print(ex)
                sys.exit(1)
            except Exception as ex:
                _log.error(ex)

from pymongo import MongoClient



for xx in orgs:
    organism = xx[0]
    db = MongoClient().pdb
    wd = "/data/organismos/" + organism + "/estructura/sndg/modelos/"
    def str_path(wd, modeldoc):
        return  glob("/".join(
            [wd, modeldoc.name + ".pdb"]))[0]
    sa = StructureAnotator(wd, struct_path= str_path)
    total = sa.total(db, organism, {})

    with tqdm(sa.iterator(db, organism, {}), total=total) as pbar:
        for model in pbar:
            pbar.set_description(model.name)

            template = model.templates[0]
            protein = Protein.objects(organism=organism, alias=template.aln_query.name).get()
            sa.annotate_model(model, protein.domains())
            model.save()

    index_seq_collection(mdb.db, organism, pathways=False, go=False, keywords=False, ec=False, organism_idx=False,
                         structure=True)


for xx in orgs:
    name = xx[0]
    load_pathways(name, "/data/organismos/" + name + "/annotation/pathways/small.sbml", mdb.db,
                  "/data/organismos/" + name + "/annotation/pathways",
                  gregexp="\(([\-\w\.]+)\)", filter_file="allfilters_con_c.dat")
    index_seq_collection(mdb.db,name,pathways=True,go=True,keywords=True,ec=True,organism_idx=True,structure=False)
    build_statistics(mdb.db,name)



import subprocess as sp
for xx in orgs:
    x = xx[0]
    for y in tqdm(glob("/data/organismos/" + x + "/estructura/sndg/modelos/*.pocket.json")):
        sp.call("ln %s %s" % (y, ("/data/organismos/%s/estructura/sndg/pockets/" % x) + y.split("/")[-1].replace(".pdb.pocket","") ),shell=True )

# rsync -avzh /data/organismos/Mbovis/estructura/sndg/pockets/* marvin:/data/organismos/Mbovis/estructura/sndg/pockets/
# rsync -avzh /data/organismos/Mbovis/estructura/sndg/modelos/*.pdb marvin:/data/organismos/Mbovis/estructura/sndg/modelos/

# ssh -L 27018:localhost:27017 marvin

pdbdb = dbdst = MongoClient(port=27017).pdb
dbdst = MongoClient(port=27018).pdb
for xx in orgs:
    x = xx[0]
    mdb.copy_genome(x,dst_db=MongoClient(port=27018).tdr)
    for s in tqdm(pdbdb.structures.find({"organism":x}),total=pdbdb.structures.count({"organism":x})):
        dbdst.structures.insert(s)
