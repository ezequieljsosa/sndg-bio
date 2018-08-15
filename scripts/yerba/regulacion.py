import networkx as nx
import pymongo
from tqdm import tqdm


def armar_grafo_completo():
    regiones = list(GFF.parse("/data/organismos/ILEX_PARA2/regulation/ncbi_IP4.gff3.reg"))
    scaffolds_dict = defaultdict(lambda: [])
    for c in regiones:
        for f in c.features:
            motif = f.sub_features[0].qualifiers["description"][0]
            if len(motifs) > 1:
                print set(motifs)
            scaffolds_dict[c.id].append({
                "start": int(f.location.start),
                "end": int(f.location.end),
                "strand": f.location.strand,
                "tf": f.qualifiers["TF_id"][0],
                "motif": motif,
                "count": len(f.sub_features)
            })
    ann = list(GFF.parse("/data/organismos/ILEX_PARA2/ncbi_IP4.gff3.whole"))
    tf_gene = []
    dist = defaultdict(lambda: {})
    contig_map = {}
    for contig in ann:
        for f in contig.features:
            if f.type == "gene" and f.sub_features:
                contig_map[f.sub_features[0].id] = contig.id

    for contig in ann:
        if contig.id in scaffolds_dict:
            for f1 in contig.features:
                for f2 in contig.features:
                    if (f1.type == "gene" and f1.sub_features) and (f2.type == "gene" and f2.sub_features):
                        dist[f1.sub_features[0].id][f2.sub_features[0].id] = abs(f1.location.start - f2.location.start)

            for reg in scaffolds_dict[contig.id]:
                if reg["strand"] == 1:
                    genes = [x for x in contig.features if int(x.location.start + 100) >= reg["end"] if
                             x.type == "gene" and x.sub_features]
                    gene = sorted(genes, key=lambda x: x.location.start)[0].sub_features[0]
                else:
                    genes = [x for x in contig.features if (x.location.end - 100) <= reg["start"] if
                             x.type == "gene" and x.sub_features]
                    try:
                        gene = sorted(genes, key=lambda x: x.location.end)[-1].sub_features[0]
                    except:
                        print([contig.id, reg["start"]])
                tf_gene.append({
                    "scaffold": contig.id, "tf": reg["tf"], "motif": reg["motif"], "count": reg["count"],
                    "gene": gene.id
                })
    nodes = {}
    edges = defaultdict(lambda: [])
    tfs = {}
    for row in tf_gene:
        nodes[row["tf"]] = 1
        tfs[row["tf"]] = [row["motif"], row["count"]]
        nodes[row["gene"]] = 1
        edges[row["tf"]].append(row["gene"])

    diG = nx.DiGraph()
    for node, scaffold in nodes.items():
        node = node.replace("*", "")
        if node == "unknown":
            continue
        diG.add_node(node, scaffold=contig_map[node],
                     count=tfs[node][1] if node in tfs else "",
                     motif=tfs[node][0] if node in tfs else "")

    for ft, genes in edges.items():
        ft = ft.replace("*", "")
        if ft == "unknown":
            continue
        for regulated in genes:
            if regulated == "unknown":
                continue
            regulated = regulated.replace("*", "")
            diG.add_edge(ft, regulated)

    assert "unknown" not in diG
    nx.write_gml(diG, "/data/organismos/ILEX_PARA2/regulation/graph.gml")


def cargar_fts_en_genoma():
    db = pymongo.MongoClient().saureus
    diG = nx.read_gml("/data/organismos/ILEX_PARA2/regulation/graph.gml")
    tfs = {x: [data["motif"], data["count"]] for x, data in diG.nodes(data=True) if "motif" in data and data["motif"]}
    print("tfs: " + str(len(tfs)))
    contig_map = {x: data["scaffold"] for x, data in diG.nodes(data=True)}
    print("contig_map: " + str(len(contig_map)))
    db.proteins.update({"organism": "ILEX_PARA2"}, {"$set": {"search.transcription_factor": "No"}}, multi=True)
    for node, data in tqdm(diG.nodes(data=True)):
        if node != "unknown":

            in_edges = diG.in_edges(node)
            out_edges = diG.out_edges(node)
            regulation = {"contig": contig_map[node], "regulated_by": []}

            if node in tfs:
                regulation["tf_type"] = data["motif"]
                regulation["regulates"] = []
                for _, reg_nofe in out_edges:
                    if reg_nofe != "unknown":
                        regulated = db.proteins.find_one({"organism": "ILEX_PARA2", "alias": reg_nofe},
                                                         {"gene": 1, "description": 1, "search.with_transcript": 1})
                        regulation["regulates"].append(
                            {"contig": contig_map[reg_nofe], "locus_tag": regulated["gene"][0],
                             "desc": regulated["description"],"is_tf": "Yes" if reg_nofe in tfs else "No",
                             "with_transcript": regulated["search"]["with_transcript"]})
            for regulates_node, _ in in_edges:
                if regulates_node != "unknown":
                    regulated = db.proteins.find_one({"organism": "ILEX_PARA2", "alias": regulates_node},
                                                     {"gene": 1, "description": 1, "search.with_transcript": 1})
                    regulation["regulated_by"].append(
                        {"contig": contig_map[regulates_node], "locus_tag": regulated["gene"][0],
                         "desc": regulated["description"],"tf_type":tfs[regulates_node][0],
                         "with_transcript": regulated["search"]["with_transcript"]})

            db.proteins.update({"organism": "ILEX_PARA2", "alias": node}, {
                "$set": {"search.transcription_factor": "Yes" if node in tfs else "No",
                         "tregulation": regulation
                         }})


cargar_fts_en_genoma()
