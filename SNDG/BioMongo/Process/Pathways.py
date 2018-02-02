def load_pathways(self, genome_name, sbml_path, prefered_biocyc=None,
                  gregexp="\(([\-\w\.]+)\)", gene_map=None, filter_file="allfilters_con_c.dat", sbml_desc=""):
    '''
    * gene_map: pikle file dict with sbml_gene_id as key and proteome_gene_id as value
    gregexp:  default = "\(([\-\w\.]+)\)"
    '''
    assert os.path.exists(sbml_path), sbml_path
    if gene_map:
        assert os.path.exists(gene_map), "%s does not exists" % gene_map
        gene_map = pickle.load(open(gene_map))


    _log.info(self.db.proteins.update({"organism":genome_name, "properties._type":"pathways"}, {"$pull":{"properties":{"_type":"pathways"}}}, multi=True))


    collection = SeqCollection.objects(name=genome_name).get()

    _log.info(self.db.proteins.update({"organism":genome_name, "properties._type":"pathways"}, {"$pull":{"properties":{"_type":"pathways"}}}, multi=True))
    _log.info(self.db.proteins.update({"organism":genome_name, "reactions":{"$exists":True}}, {"$set":{"reactions":[]}}, multi=True))


    ps = PathwaysAnnotator(self.db, genome_name, self.fs_resolver.pathways_dir(genome_name))
    ps.sbml(sbml_path)
    ps.species_filter(filter_file)
    gene_name_regexp = re.compile(gregexp)
    if gene_map:
        def map_gene(notes):
            genes = gene_name_regexp.findall(notes)
            return [gene_map[x] if x in gene_map else x for x in genes]
        ps.extract_genes_from_notes(map_gene)
    else:
        ps.extract_genes_from_notes(lambda notes : gene_name_regexp.findall(notes))
    if prefered_biocyc:
        ps.prefered_biocyc = prefered_biocyc
    ps.annotate()

    result = AnnotationPipelineResult(name="pathways", version=1, timestamp=datetime.datetime.now)
    result.inputs = {
        "prefered_biocyc":prefered_biocyc,
        "gregexp":gregexp,
        "gene_map":gene_map,
        "sbml":gene_map,
        "sbml_desc":sbml_desc
    }

    chokepoints = self.db.proteins.count({"organism":genome_name, "properties.property":"chokepoint"})

    result.results = {
        "pathways": len(ps.sbmlprocessor.pathways) ,
        "annotated_proteins": self.db.proteins.count({"organism":genome_name, "reactions.0":{"$exists":True}}) ,
        "chokepoints":chokepoints,
        "unmapped_genes":ps.unmapped_genes
    }

    _log.info("PW annotation finished")

    collection.pathways = ps.pathways
    collection.pipelines.append(result)

    collection.save()

def correct_chokes(self, name):

    metabolites_in = defaultdict(list)
    metabolites_out = defaultdict(list)
    for p in self.db.proteins.find({"organism":name, "reactions.0" :{"$exists": True}}):
        for r in p["reactions"]:
            for m in r["products"]:
                metabolites_out[m["name"]].append(r["name"])
            for m in r["substrates"]:
                metabolites_in[m["name"]].append(r["name"])

    for m, r in metabolites_in.items():
        if (len(set(r)) > 1):
            # or (self.db.proteins.count({"organism":name,"reactions.name": r[0]}) > 1)):
            del metabolites_in[m]
    for m, r in metabolites_out.items():
        if (len(set(r)) > 1):
            # or (self.db.proteins.count({"organism":name,"reactions.name": r[0]}) > 1)):
            del metabolites_out[m]

    choke_reactions_in = []
    for rs in metabolites_in.values():
        choke_reactions_in += rs
    choke_reactions_out = []
    for rs in metabolites_out.values():
        choke_reactions_out += rs

    reaction_metabolites = defaultdict(lambda : [])
    for m, rs in metabolites_in.items():
        for r in rs:
            reaction_metabolites[r].append(m)
    for m, rs in metabolites_out.items():
        for r in rs:
            reaction_metabolites[r].append(m)

    for p in Protein.objects(organism=name, reactions__0__exists=True).no_cache().timeout(False):
        cout = bool([r.name for r in p.reactions if r.name in  choke_reactions_out])
        cin = bool([r.name for r in p.reactions if r.name in  choke_reactions_in])
        p.search.chokepoint = cout | cin
        if p.search.chokepoint:
            p.search.chokepoint_type = "double" if (cout & cin) else ("production" if cout else "consuming")
            prop = [x for x in p.properties if x.property == "chokepoint"]
            if prop:
                prop = prop[0]
                prop.metabolites = []
                for x in p.reactions:
                    if x.name in reaction_metabolites:
                        prop.metabolites += reaction_metabolites[x.name]
            else:
                prop = BioProperty(_type="pathways", property="chokepoint", metabolites=[] , type=p.search.chokepoint_type)
                for x in p.reactions:
                    if x.name in reaction_metabolites:
                        prop.metabolites += reaction_metabolites[x.name]
                p.properties.append(prop)
        else:
            del p.search.chokepoint_type
            p.properties = [x for x in p.properties if x.property != "chokepoint"]
        p.save()