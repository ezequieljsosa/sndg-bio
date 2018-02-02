def annotate_variants(self, organism_name, strain_name, database, parse_change):
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


def annotate_variants_with_prots(self, organism_name, dbs, drugs, force=False):
    """
    drugs: list of strings, example TBDream.drugs or Saureus.drugs
    """
    for idx, p in enumerate(Protein.objects(
            __raw__={"organism": organism_name, "features.qualifiers.strain": {"$exists": 1}}).no_cache()):
        print idx
        pvariants = list(VarDoc.objects(organism=organism_name, gene__in=p.gene))

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


def indexVariants(self, organism):
    print self.db.var_col_ont_idx.remove({"seq_collection_name": organism})
    for ont in SeqColOntologyIndex.objects(seq_collection_name=organism):
        ont.id = ObjectId()
        ont.count = VarDoc.objects(organism=organism, ontologies=ont.term).count()
        vidx = VarSeqColOntologyIndex(**ont._data)
        if ont.count:
            vidx.save()
