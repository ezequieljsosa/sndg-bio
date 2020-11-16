from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
from SNDG.BioMongo.Model.SeqCollection import SeqCollection, SeqColDruggabilityParam
from SNDG.WebServices.Offtargeting import Offtargeting
from SNDG import init_log, mkdir, execute
from SNDG.WebServices import PROXIES
import os

PROXIES["ftp_proxy"] = "http://proxy.fcen.uba.ar:8080"
init_log()

mdb = BioMongoDB("tdr",port=27018)

off_props = {"human_offtarget": SeqColDruggabilityParam(**{
    "target": "protein",
    "defaultGroupOperation": "max",
    "defaultValue": 0.6,
    "name": "human_offtarget",
    "defaultOperation": ">",
    "_cls": "SeqColDruggabilityParam",
    "uploader": "demo",
    "_class": "ar.com.bia.entity.SeqCollectionDoc",
    "type": "number",
    "options": [],
    "description": "This score reflects the results of a blastp search of the pathogen protein in the human proteome database (Gencode v17) with the scale 1 - max(alignment identity), so when a protein has no hit in the human proteome, the value is 1, and if it has 2 hits, one with an identity of 0.4 and other with 0.6, the score is 0.4 (human_offtarget = 1 - 0.6, uses the max identity)."
}),
             "gut_microbiome": SeqColDruggabilityParam(**{
                 "target": "protein",
                 "name": "gut_microbiome",
                 "_cls": "SeqColDruggabilityParam",
                 "uploader": "demo",
                 "_class": "ar.com.bia.entity.SeqCollectionDoc",
                 "type": "number",
                 "options": [],
                 "description": "Number of gut microbiome organisms that have at least one hit (blast: identity > 40% evalue 1e-5)"
             }),
             "gut_microbiome_norm": SeqColDruggabilityParam(**{
                 "target": "protein",
                 "name": "gut_microbiome_norm",
                 "_cls": "SeqColDruggabilityParam",
                 "uploader": "demo",
                 "_class": "ar.com.bia.entity.SeqCollectionDoc",
                 "type": "number",
                 "options": [],
                 "description": "gut_microbiome normalized by the total number of compared bacteria (226 - https://doi.org/10.1038/s41598-018-28916-7 - Supplementary Table S1)"
             }),
             "hit_in_deg": SeqColDruggabilityParam(**{
                 "target": "protein",
                 "defaultGroupOperation": "avg",
                 "defaultValue": "Yes",
                 "name": "hit_in_deg",
                 "defaultOperation": "equal",
                 "_cls": "SeqColDruggabilityParam",
                 "uploader": "demo",
                 "_class": "ar.com.bia.entity.SeqCollectionDoc",
                 "type": "value",
                 "options": [
                     "No",
                     "Yes"
                 ],
                 "description": "Has a hit in Database of Essential Genes"
             })

             }
from SNDG.Sequence import read_blast_table
from tqdm import tqdm

# cols = list(SeqCollection.objects(name__nin=["cruzi","pdb"]))
cols = list(SeqCollection.objects(name__nin=["cruzi","pdb"]))
cpus = 4
db = mdb.db
for seqCol in tqdm(cols):
    mkdir("/data/organismos/" + seqCol.name + "/contigs")
    proteome = "/data/organismos/" + seqCol.name + "/contigs/genoma.fasta"
    if not os.path.exists(proteome):
        mdb.protein_fasta(proteome, seqCol.name)

    out = "/data/organismos/" + seqCol.name + "/annotation/offtarget/"
    mkdir(out)

    if not seqCol.has_druggability_param("human_offtarget"):

        seqCol.druggabilityParams.append(off_props["human_offtarget"])
        db = "/data/databases/human/gencode.v17.pc_translations.fa"

        execute(
            "blastp -evalue 1e-5 -max_hsps 1 -outfmt 6 -max_target_seqs 1 -db {db} -query {query} -out {out} -num_threads {cpus}",
            db=db, query=proteome, out=out + "human_offtarget.tbl", cpus=cpus)

        mdb.db.proteins.update({"organism": seqCol.name},
                               {"$set": {"search.human_offtarget": 1}}, multi=True)
        for _, r in tqdm(read_blast_table(out + "human_offtarget.tbl").iterrows()):
            mdb.db.proteins.update({"organism": seqCol.name, "gene": r["query"]},
                                   {"$set": {"search.human_offtarget": 1 - r.identity / 100}})



    if not seqCol.has_druggability_param("gut_microbiome"):
        seqCol.druggabilityParams.append(off_props["gut_microbiome"])
        seqCol.druggabilityParams.append(off_props["gut_microbiome_norm"])
        db = "/data/databases/human/gut_microbiota.fasta"

        execute(
            "blastp -evalue 1e-5 -max_hsps 1 -outfmt 6  -db {db} -query {query} -out {out} -num_threads {cpus}",
            db=db, query=proteome, out=out + "gut_microbiome.tbl", cpus=cpus)

        mdb.db.proteins.update({"organism": seqCol.name},
                               {"$set": {"search.gut_microbiome": 0,
                                         "search.gut_microbiome_norm": 0
                                         }}, multi=True)

        prot_off = Offtargeting.count_organism_from_microbiome_blast(out + "gut_microbiome.tbl", db)
        for locus_tag, organisms in tqdm(prot_off.items()):
            assert len(organisms) <= 226
            mdb.db.proteins.update({"organism": seqCol.name, "gene": locus_tag},
                                   {"$set": {"search.gut_microbiome": len(organisms),
                                             "search.gut_microbiome_norm": len(organisms) / 226.0
                                             }})



    if not seqCol.has_druggability_param("hit_in_deg"):
        seqCol.druggabilityParams.append(off_props["hit_in_deg"])

        db = "/data/databases/deg/degaa-p.dat"
        if seqCol.name in ["cruzi", "LMajor", "SC_MAS", "PVIVAX", "TGONDII"]:
            db = "/data/databases/deg/degaa-e.dat"

        execute(
            "blastp -evalue 1e-5 -max_hsps 1 -qcov_hsp_perc 80 -outfmt 6 -db {db} -query {query} -out {out} -num_threads {cpus}",
            db=db, query=proteome, out=out + "deg.tbl", cpus=cpus)

        mdb.db.proteins.update({"organism": seqCol.name},
                               {"$set": {"search.hit_in_deg": "No"}}, multi=True)
        for _, r in tqdm(read_blast_table(out + "deg.tbl").iterrows()):
            if (r.identity / 100.0) > 0.7:
                mdb.db.proteins.update({"organism": seqCol.name, "gene": r["query"]},
                                       {"$set": {"search.hit_in_deg": "Yes"}})

    seqCol.save()
    from SNDG.BioMongo.Process.Index import build_statistics
    build_statistics(mdb.db,seqCol.name)