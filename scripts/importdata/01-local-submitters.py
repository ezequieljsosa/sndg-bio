def obtener_proyectos_argentinos(self):
    all_arg_proj = Entrez.read(Entrez.esearch(db="bioproject", term="argentina", retmax=10000))
    _log.info("existen %s proyectos con la palabra argentina" % all_arg_proj["Count"])
    for pid in tqdm(all_arg_proj["IdList"]):
        yield Entrez.read(Entrez.esummary(db="bioproject", id=pid))["DocumentSummarySet"]["DocumentSummary"][0]

def obtener_submitters(self):
    submitters = []
    for bioproject in self.obtener_proyectos_argentinos():
        submitters = submitters + bioproject["Submitter_Organization_List"]

    return set(submitters)

def actualizar_submitters(self, submitters):

    for submitter_raw in tqdm(submitters):
        submitter = str(submitter_raw.encode("utf-8"))
        submitter_in_db = Submitter.select().where((Submitter.source == "ncbi") & (Submitter.name == submitter))
        if not submitter_in_db.exists():
            submitter_in_db = Submitter(name=submitter, source="ncbi")
            submitter_in_db.save()

def update_assemblies(self, submitters):
    db = "assembly"
    for submitter in tqdm(submitters):
        data = Entrez.read(Entrez.esearch(db=db, term='"' + submitter.name + '"[Submitter Organization]'))
        for recurso_id in data["IdList"]:
            er_in_db = ExternalAssembly.select().where(ExternalAssembly.identifier == recurso_id)
            if not er_in_db.exists():
                recurso = Entrez.read(Entrez.esummary(db=db, id=recurso_id))
                data = self.resource_handler[db].attributes(recurso)
                name = str(self.resource_handler[db].name(recurso))
                genome = str(data["SpeciesName"])
                ea = ExternalAssembly(type=db, name=name, identifier=recurso_id
                                      , assembly_accession=data['AssemblyAccession'], genome=genome
                                      , assembly_name=data['AssemblyName'])
                ea.save(force_insert=True)
            else:
                ea = list(er_in_db)[0]

            AssemblySubmitters.create(submitter=submitter, resource=ea)

def actualizar_recursos_externos(self):
    submitters = Submitter.select().where((Submitter.source == "ncbi") & (Submitter.rejected == False))
    for submitter in submitters:
        for db in NCBI.dbs_con_submitter:
            data = Entrez.read(Entrez.esearch(db=db, term='"' + submitter.name + '"[Submitter Organization]'))
            if 'ErrorList' in data:
                _log.warn("error buscando %s en %s : %s " % (submitter.name, db, str(data['ErrorList'])))
            if db == "bioproject":
                bioprojects = data["IdList"]
            if db == "assembly":
                assemblies = data["IdList"]
            for recurso_id in data["IdList"]:
                er_in_db = ExternalResource.select().where((ExternalResource.submitter == submitter) &
                                                           (ExternalResource.identifier == recurso_id))
                if not er_in_db.exists():
                    recurso = Entrez.read(Entrez.esummary(db=db, id=recurso_id))
                    nombre = str(self.resource_handler[db].name(recurso))
                    ExternalResource(submitter=submitter, type=db, name=nombre, identifier=recurso_id).save()
                    if db == "assembly":
                        ExternalAssembly(submitter=submitter, type=db, name=nombre, identifier=recurso_id
                                         , assembly_accession=recurso['Assembly_Accession']
                                         , assembly_name=recurso['Assembly_Name']).save()

        for db in (set(NCBI.dbs) - set(NCBI.dbs_con_submitter)):
            for bioproject in bioprojects:
                data = Entrez.read(Entrez.esearch(db=db, term=bioproject + '[BioProject]'))
                if 'ErrorList' in data and not data['ErrorList']['PhraseNotFound']:
                    _log.warn("error buscando %s en %s : %s " % (submitter, db, str(data['ErrorList'])))
                for recurso_id in data["IdList"]:
                    er_in_db = ExternalResource.select().where((ExternalResource.submitter == submitter) &
                                                               (ExternalResource.identifier == recurso_id))
                    if not er_in_db.exists():
                        recurso = Entrez.read(Entrez.esummary(db=db, id=recurso_id))
                        nombre = self.resource_handler[db].name(recurso).encode("utf-8")
                        ExternalResource(submitter=submitter, type=db, name=nombre, identifier=recurso_id).save()

        db = "genome"
        for assembly in assemblies:
            data = Entrez.read(Entrez.esearch(db=db, term=assembly + '[AssemblyID]'))
            if 'ErrorList' in data and not data['ErrorList']['PhraseNotFound']:
                _log.warn("error buscando %s en %s : %s " % (submitter, db, str(data['ErrorList'])))
            for recurso_id in data["IdList"]:
                er_in_db = ExternalResource.select().where((ExternalResource.submitter == submitter) &
                                                           (ExternalResource.identifier == recurso_id))
                if not er_in_db.exists():
                    recurso = Entrez.read(Entrez.esummary(db=db, id=recurso_id))
                    nombre = self.resource_handler[db].name(recurso).encode("utf-8")
                    # 'Assembly_Name': 'v.1.0'
                    # 'Assembly_Accession': 'GCA_000004255.1'
                    ExternalResource(submitter=submitter, type=db, name=nombre, identifier=recurso_id).save()