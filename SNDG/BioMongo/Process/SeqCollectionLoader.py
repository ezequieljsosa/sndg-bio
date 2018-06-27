"""

"""


class SeqCollectionLoader:
    def load_genome_fasta_and_gff(self, genome_name, fasta_path, gff_path, organism="", description=""):
        from BIADeploy.BioDocFactory import BioDocFactory
        g = Genome(
            name=genome_name,
            description=description,
            organism=organism,
            auth=str(BioMongoDB.demo_id))
        g.save()

        #         gene_ids = {}
        contig_dict = {}
        for i, seqrecord in enumerate(bpio.parse(fasta_path, "fasta")):
            _log.debug("contig %s %i" % (seqrecord.id, i))

            contig_dict[seqrecord.id] = seqrecord
        #             gene_ids.update(gene_ids2)

        with open(gff_path) as h:
            gffs = h.read().split("##gff-version 3")[1:]
        for gff in gffs:
            records = list(GFF.parse(StringIO.StringIO(gff)))
            prots = []
            for i, r in enumerate(records):
                if r.features:
                    r.seq = contig_dict[r.id].seq
                    contig, _ = BioDocFactory.create_contig(r, g, source="SNDG_pipeline")  # gene_ids2
                    contig.save()
                for f1 in r.features:

                    if f1.type == "gene":
                        for f in f1.sub_features:
                            if f.type == "mRNA":
                                protDoc = BioDocFactory.create_protein(contig_dict[r.id], f,
                                                                       [e for e in f.sub_features if e.type == "exon"])
                                if len(protDoc.seq) > 30000:
                                    raise Exception("No existen proteinas tan largas...")
                                # protDoc.gene_id = gene_ids[f.qualifiers["locus_tag"][0]]
                                protDoc.organism = genome_name
                                protDoc.auth = str(BioMongoDB.demo_id)
                                protDoc.seq_collection_id = g
                                prots.append(protDoc)
                                if i and ((i % 1000) == 0):
                                    print "se procesaron %i proteinas" % i
                                    Protein.objects.insert(prots)
                                    prots = []
                if prots:
                    Protein.objects.insert(prots)
                    prots = []

    def explore_gff_from_genome(self, gname):
        return self.explore_gff(self.fs_resolver.gff_annotation(gname))

    def explore_gff(self, gff_path):

        from BCBio.GFF import GFFExaminer
        examiner = GFFExaminer()
        with open(gff_path) as h:
            parentchild = examiner.parent_child_map(h)
            pprint.pprint(parentchild)
        with open(gff_path) as h:
            pprint.pprint(examiner.available_limits(h))

    #     GFF.parse("./Staphylococcus_aureus_subsp_aureus_n315.ASM964v1.35.gff3")

    def travel_gff(self, records):
        data = defaultdict(lambda: [])
        examples = {}

        def travel_feature(feature, ftype):
            if feature.sub_features:
                for f in feature.sub_features:
                    path = ftype + "|" + f.type
                    travel_feature(f, path)
                    data[path] = list(set(f.qualifiers.keys() + data[path]))
                    if path in examples:
                        examples[path].update(f.qualifiers)
                    else:
                        examples[path] = f.qualifiers

        for r in records:
            for f in r.features:
                data[f.type] = list(set(f.qualifiers.keys() + data[f.type]))
                travel_feature(f, f.type)

        for x, y in data.items():
            print x
            if x in examples:
                for h in y:
                    print "->" + h + ": " + str(examples[x][h])
            else:
                print y

    def genome_from_fasta(self, name, organism, dataDir, fasta, auth,
                          swiss_prot="./annotation/uniprot/sp.xml",
                          trembl="./annotation/uniprot/tr.xml",
                          swiss_prot_blast="./annotation/uniprot/sp_blast.xml",
                          trembl_blast="./annotation/uniprot/tr_blast.xml"):

        os.chdir(dataDir)
        self.delete_seq_collection(name)
        genome = Genome(name=name, organism=organism, type="genome", auth=str(auth))
        genome.save()
        for seq in bpio.parse(fasta, "fasta"):
            Protein(seq=str(seq.seq), gene=list([seq.id]), name=seq.id,
                    description=seq.description, ontologies=[],
                    features=[], organism=name,
                    seq_collection_id=genome.id,
                    alias=list(set([seq.id] + [seq.id.split("|")[0]] if "|" in seq.id else []))).save()

        uannotator = UniprotAnnotator()
        uannotator.auth(auth)
        uannotator.name(name)
        uannotator.fasta(fasta)
        uannotator.organism(organism)

        uannotator.swiss_prot(swiss_prot)
        uannotator.trembl(trembl)

        uannotator.swiss_prot_blast(swiss_prot_blast)
        uannotator.trembl_blast(trembl_blast)

        uannotator.annotate_proteome(Protein.objects(organism=name).no_cache().timeout(False),
                                     uannotator.get_homologs())
