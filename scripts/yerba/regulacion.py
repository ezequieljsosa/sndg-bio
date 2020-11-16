import networkx as nx
import pymongo
from tqdm import tqdm
from BCBio import GFF
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation


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
    # dist = defaultdict(lambda: {})
    contig_map = {}
    for contig in ann:
        for f in contig.features:
            if f.type == "gene" and f.sub_features:
                contig_map[f.sub_features[0].id] = contig.id

    for contig in ann:
        if contig.id in scaffolds_dict:
            # for f1 in contig.features:
            #     for f2 in contig.features:
            #         if (f1.type == "gene" and f1.sub_features) and (f2.type == "gene" and f2.sub_features):
            #             dist[f1.sub_features[0].id][f2.sub_features[0].id] = abs(f1.location.start - f2.location.start)

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
    diGString = nx.read_gml("/data/organismos/ILEX_PARA2/regulation/graph_string60.gml").to_undirected()

    tfs = {x: [data["motif"], data["count"]] for x, data in diG.nodes(data=True) if "motif" in data and data["motif"]}
    print("tfs: " + str(len(tfs)))
    contig_map = {x: data["scaffold"] for x, data in diG.nodes(data=True)}
    print("contig_map: " + str(len(contig_map)))
    db.proteins.update({"organism": "ILEX_PARA2"},
                       {"$set": {"search.transcription_factor": "No"}, "$unset": {"tregulation": ""}}, multi=True)
    for node, data in tqdm(diG.nodes(data=True)):
        if node != "unknown":

            in_edges = diG.in_edges(node)
            out_edges = diG.out_edges(node)
            regulation = {"contig": contig_map[node], "regulated_by": []}
            if node in diG:
                if node in tfs:
                    regulation["tf_type"] = data["motif"]
                if node in diGString:
                    regulation["string"] = diGString.node[node]["string"]

                regulation["regulates"] = []
                for _, reg_nofe in out_edges:
                    if reg_nofe != "unknown":
                        regulated = db.proteins.find_one({"organism": "ILEX_PARA2", "alias": reg_nofe},
                                                         {"gene": 1, "description": 1, "search.with_transcript": 1})

                        score = (diGString.edges[node, reg_nofe]["score"][0]
                        if diGString.has_edge(node, reg_nofe) else 0)
                        if score:
                            if "string" not in regulation:
                                raise Exception()
                        regulation["regulates"].append(
                            {"contig": contig_map[reg_nofe], "locus_tag": regulated["gene"][0],
                             "string": diGString.node[reg_nofe]["string"] if reg_nofe in diGString else None,
                             "score": score,
                             "desc": regulated[
                                 "description"], "is_tf": "Yes" if reg_nofe in tfs else "No",
                             "with_transcript":
                                 regulated["search"]["with_transcript"]})
                for regulates_node, _ in in_edges:
                    if regulates_node != "unknown":
                        regulated = db.proteins.find_one({"organism": "ILEX_PARA2", "alias": regulates_node},
                                                         {"gene": 1, "description": 1, "search.with_transcript": 1})

                        score = (diGString.edges[regulates_node, node]["score"][0]
                        if diGString.has_edge(node, regulates_node) else 0)

                        if score:
                            if "string" not in regulation:
                                raise Exception()

                        regulation["regulated_by"].append(
                            {"contig": contig_map[regulates_node], "locus_tag": regulated["gene"][0],
                             "string": diGString.node[regulates_node][
                                 "string"] if regulates_node in diGString else None,
                             "score": score,
                             "desc": regulated["description"], "tf_type": tfs[regulates_node][0],
                             "with_transcript": regulated["search"]["with_transcript"]})

                db.proteins.update({"organism": "ILEX_PARA2", "alias": node}, {
                    "$set": {"search.transcription_factor": "Yes" if node in tfs else "No",
                             "tregulation": regulation
                             }})


def regulado_por():
    db = pymongo.MongoClient().saureus
    diG = nx.read_gml("/data/organismos/ILEX_PARA2/regulation/graph.gml")

    dp = {
        "target": "protein",
        "defaultGroupOperation": "avg",
        "defaultValue": "Yes",
        "name": "regulated_by_tf",
        "defaultOperation": "equal",
        "_cls": "SeqColDruggabilityParam",
        "uploader": "demo",
        "_class": "ar.com.bia.entity.SeqCollectionDoc",
        "type": "value",
        "options": [
            "Yes",
            "No"
        ],
        "description": "Gene is potentially regulated by a known TF"
    }
    print (db.sequence_collection.update({"name": "ILEX_PARA"}, {"$addToSet": {"druggabilityParams": dp}}))
    print(db.proteins.update({"organism": "ILEX_PARA"}, {"$set": {"search.regulated_by_tf": "No"}}))

    for node in tqdm(diG.nodes()):
        if diG.in_degree(node):
            db.proteins.update({"organism": "ILEX_PARA", "alias": node},
                               {"$set": {"search.regulated_by_tf": "Yes"}}
                               )


def agrupar_sitios():
    regiones = list(GFF.parse("/data/organismos/ILEX_PARA2/regulation/ncbi_IP4.gff3.reg"))
    ids = 1
    groups = {}
    for c in tqdm(regiones):
        groups[c.id] = []
        for strand in [1,-1]:
            group = SeqFeature(id=c.features[0], type="grouped_transcription_regulatory_region",
                               location=c.features[0].location)
            group.sub_features = []

            fs = sorted([f for f in c.features if f.strand == strand],key=lambda x:x.location.start)
            if not fs:
                continue
            group.sub_features += [fs[0]]

            for f in fs[1:]:
                end = max([x.location.end for x in group.sub_features])
                if (( abs(f.location.start - end) < 1500)
                     or (set(range(f.location.start, f.location.end)) &
                         set(range(group.sub_features[-1].location.start, group.sub_features[-1].location.end)))
                ) :
                    group.sub_features.append(f)
                else:
                    group.qualifiers = {
                        "description": "_".join(sorted(set([x.qualifiers["description"][0].split(" regulatory region")[0]
                                                 for x in group.sub_features]))), "ID": ["ILEXPARARR" + str(ids)]}
                    ids += 1

                    group.location = FeatureLocation(start=min([x.location.start for x in group.sub_features]),
                                                     end=max([x.location.end for x in group.sub_features]),
                                                     strand=f.location.strand)
                    assert group.location.start < group.location.end
                    if (group.location.end -group.location.start) > 5000:
                        print(group.qualifiers["ID"])


                    groups[c.id].append(group)
                    group = SeqFeature(id=c.features[0], type="grouped_transcription_regulatory_region",
                                       location=f.location)
                    group.sub_features = [f]
            if group:
                group.qualifiers = {"description": "_".join(sorted(set([x.qualifiers["description"][0].split(" binding site")[0]
                                                             for x in group.sub_features]))),
                                    "ID": ["ILEXPARARR" + str(ids)]}
                ids += 1
                group.location = FeatureLocation(start=min([x.location.start for x in group.sub_features]),
                                                 end=max([x.location.end for x in group.sub_features]),
                                                 strand=f.location.strand)
                assert group.location.start < group.location.end
                if (group.location.end -group.location.start) > 5000:
                    print(group.qualifiers["ID"])

                groups[c.id].append(group)


    # for _, v in groups.items():
    #     for x in v:
    #         x.sub_features = []
    records = [SeqRecord(id=k, name="", description="", seq=Seq(""), features=v)
                 for k, v in groups.items()]
    GFF.write(tqdm(records), open("/data/organismos/ILEX_PARA2/regulation/grouped.gff", "w"))


agrupar_sitios()
# cargar_fts_en_genoma()


# PERL5LIB=/home/eze/perl5/lib/perl5 time ./bin/flatfile-to-json.pl --gff "/data/organismos/ILEX_PARA2/regulation/ncbi_IP4.gff3.reg" --out "./data/ILEX_PARA2" --key "Regulation Sites" --trackLabel "Regulation" --trackType FeatureTrack

# {
#     "compress" : 0,
#     "key" : "Regulation Sites",
#     "label" : "Regulation",
#     "storeClass" : "JBrowse/Store/SeqFeature/NCList",
#     "style" : {
#         "color": "red",
#         "label":""
#     },
#     "trackType" : "FeatureTrack",
#     "type" : "FeatureTrack",
#     "urlTemplate" : "tracks/Regulation sites/{refseq}/trackData.json"
# }
