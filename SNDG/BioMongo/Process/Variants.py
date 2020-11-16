from SNDG.BioMongo.Model.Variant import Variant, Allele, SampleAllele
from SNDG.BioMongo.Model.SeqColOntologyIndex import VarSeqColOntologyIndex

from SNDG.BioMongo.Model.Protein import Protein
from SNDG.BioMongo.Model.SeqColDruggabilityParam import SeqColDruggabilityParam, SeqColDruggabilityParamTypes


def annotate_variants(organism_name, strain_name, database, parse_change):
    """
    parse_change: function that transforms  dbvar.qualifiers["change"] into aa_ref, aa_alt
    """

    collection = SeqCollection.objects(name=organism_name).get()

    prop = strain_name + "_" + database
    Protein.objects(organism=organism_name).update(__raw__={"$set": {"search." + prop: False}})
    if not collection.has_druggability_param(prop):
        dp = SeqColDruggabilityParam(name=prop,
                                     description="Variant in strain " + strain_name + " is reported in " + database,
                                     target="variant-strain",
                                     type=SeqColDruggabilityParamTypes.value, uploader="demo")
        dp.options = ["true", "false"]
        dp.defaultValue = "true"
        dp.defaultOperation = "equal"
        dp.defaultGroupOperation = "avg"
        collection.druggabilityParams.append(dp)
    prop = strain_name + "_" + database + "_pos"
    Protein.objects(organism=organism_name).update(__raw__={"$set": {"search." + prop: False}})
    if not collection.has_druggability_param(prop):
        dp = SeqColDruggabilityParam(name=prop,
                                     description="The position of the variant the strain " + strain_name + " is reported in " + database,
                                     target="variant-strain", type=SeqColDruggabilityParamTypes.value, uploader="demo")
        dp.options = ["true", "false"]
        dp.defaultValue = "true"
        dp.defaultOperation = "equal"
        dp.defaultGroupOperation = "avg"
        collection.druggabilityParams.append(dp)
    collection.save()

    for p in Protein.objects(__raw__={"organism": organism_name, "features.qualifiers.strain": strain_name}).no_cache():
        dbvars = [f for f in p.features if f.type == database]
        dirty = False
        if dbvars:
            strainvars = [f for f in p.features if
                          (f.type == "strain_variant") and (f._data["qualifiers"]["strain"] == strain_name)]

            for dbvar in dbvars:
                dirty = True
                for strainvar in strainvars:
                    strainvar._data["qualifiers"]["ref_pos"] = False
                    if dbvar.location.start == strainvar.location.start:
                        p.search[strain_name + "_" + database + "_pos"] = True
                        strainvar._data["qualifiers"]["ref_pos"] = dbvar._id

                        try:
                            dref, dalt = parse_change(dbvar._data["qualifiers"]["change"])
                            sref, salt = strainvar._data["qualifiers"]["change"].strip().split("/")
                            sref = sref.strip()
                            salt = salt.strip()
                            if (dref == sref) and (dalt == salt):
                                p.search[strain_name + "_" + database] = True
                                strainvar._data["qualifiers"]["ref"] = dbvar._id

                        except Exception as ex:
                            _log.warn(ex)
                        if (("frameshift" in dbvar._data["qualifiers"]["change"].lower())
                                and ("frameshift" in strainvar._data["qualifiers"]["change"].lower())
                        ):
                            p.search[strain_name + "_" + database] = True
            if dirty:
                p.save()


def annotate_variants_with_prots(organism_name, dbs, drugs, force=False):
    """
    drugs: list of strings, example TBDream.drugs or Saureus.drugs
    """
    for idx, p in enumerate(Protein.objects(
            __raw__={"organism": organism_name, "features.qualifiers.strain": {"$exists": 1}}).no_cache()):
        print idx
        pvariants = list(Variant.objects(organism=organism_name, gene__in=p.gene))

        for vd in pvariants:

            if ((vd.search == None) or force):

                psearch = p.search
                del psearch.structures
                vd.search = psearch
                vd.ontologies = p.ontologies

                for r in p.reactions:
                    for pw in r.pathways:
                        vd.ontologies.append(pw)
                vd.ontologies = list(set(vd.ontologies))

            for drug in drugs:
                vd.search[drug] = False
            vd.search["resistance"] = False
            for db in dbs:

                for sample_allele in vd.sample_alleles:
                    aa_pos = sample_allele.aa_pos
                    feature = [f for f in p.features if (f.type == db)
                               and (f.location.start == aa_pos)]
                    if feature:
                        feature = feature[0]
                        if feature._data["qualifiers"]:
                            vd.search["resistance"] = True

                        feature = [f for f in p.features if (f.type == db)
                                   and (f.location.start == aa_pos) and ("mut" in f._data["qualifiers"])
                                   and
                                   ((f._data["qualifiers"]["mut"] == sample_allele.aa_alt)
                                    or
                                    ((f._data["qualifiers"]["change"].lower() == "frameshift"
                                      and sample_allele.variant_type == "frameshift_variant"
                                      ))
                                    )
                                   ]
                        if feature:
                            sample_allele.feature = feature[0]
                            if sample_allele.feature._data["qualifiers"]:
                                vd.search[sample_allele.feature._data["qualifiers"]["drug"]] = True  #
            vd.save()


def indexVariants(db, organism):
    print db.var_col_ont_idx.remove({"seq_collection_name": organism})
    for ont in SeqColOntologyIndex.objects(seq_collection_name=organism):
        ont.id = ObjectId()
        ont.count = Variant.objects(organism=organism, ontologies=ont.term).count()
        vidx = VarSeqColOntologyIndex(**ont._data)
        if ont.count:
            vidx.save()


if __name__ == '__main__':
    import pymongo
    import os
    from SNDG.Comparative.VcfSnpeffIO import VcfSnpeffIO
    from SNDG.BioMongo.Process.BioMongoDB import BioMongoDB
    from tqdm import tqdm
    import subprocess
    from bson import ObjectId

    mdb = BioMongoDB("tdr")

    base = "/data/projects/23staphylo/processed/core-mapping-n315/"
    ref = "SaureusN315"

    # ref = "H37Rv"
    # base = "/data/projects/mtbxdr/processed/h37rv_processing/"

    for strain in tqdm(os.listdir(base)):
        #vcf_file0 = base + strain + "/raw_snps.vcf"
        vcf_file0 = base + strain + "/variants.vcf"
        vcf_file = base + strain + "/variants.ann.vcf"
        if os.path.exists(vcf_file0):
            # if not os.path.exists(vcf_file):
            #     subprocess.check_output("java -jar /opt/snpEff/snpEff.jar h37rv " + vcf_file0 + " > " + vcf_file   ,shell=True)

            total = int(subprocess.check_output("grep -v ^#  " + vcf_file + " | wc -l", shell=True).split()[0])
            for v, effects in tqdm(VcfSnpeffIO.parse(vcf_file), total=total):
                effect = effects[0]
                variant = Variant.objects(organism=ref, pos=v.POS, contig=v.CHROM,
                                          ref=v.REF).first()
                sample_allele = SampleAllele(sample=strain)
                allele = Allele(_id=ObjectId(), samples=[sample_allele], alt=str(v.ALT[0]),
                                variant_type=list(effect.annotation))

                for k, var_value in v.samples[0].data._asdict().items():
                    val = "|".join(map(str, var_value)) if isinstance(var_value, (list, tuple)) else str(var_value)
                    sample_allele.annotations[k] = val
                sample_allele.annotations["qual"] = v.QUAL

                if not variant:
                    variant = Variant(organism=ref, pos=v.POS, contig=v.CHROM,
                                      ref=v.REF, sample_alleles=[allele])
                    variant.search.core = True

                elif not [x for x in variant.sample_alleles if strain in [y.sample for y in x.samples]]:
                    alleledb = [x for x in variant.sample_alleles if x.alt == allele.alt]
                    if alleledb:
                        alleledb[0].samples.append(sample_allele)
                    else:
                        variant.sample_alleles.append(allele)

                else:
                    # sample allele is present
                    continue
                variant.gene = effect.geneid


                p = mdb.db.proteins.find_one({"organism":ref, "gene":effect.geneid},{"gene":1})
                if p:
                    variant.prot_ref = p["_id"]
                    variant.gene = p["gene"][0]
                else:

                    print "%s does not exists in the db" % effect.geneid

                allele.aa_ref = effect.aa_ref
                allele.aa_alt = effect.aa_alt
                allele.aa_pos = effect.aa_pos


                variant.save()
