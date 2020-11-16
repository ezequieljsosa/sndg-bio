"""
#antismash = {}
    ...: #with open("")
    ...: #for x in bpio.parse("antismash.gbk","gb"):
    ...: #    for f in x.features:
    ...: #        if "protein_id" in f.qualifiers:
    ...: #            antismash[f.qualifiers["protein_id"][0].replace("ncbi_","")] = f.qualifiers
    ...: cols = "query_name^Iseed_eggNOG_ortholog^Iseed_ortholog_evalue^Iseed_ortholog_score^Ibest_tax_level^IPreferred_name^IGOs^IEC^IKEGG_ko^IKEGG_Pathway^IKEGG_Module^IKEGG_Reaction^IKEGG_rclass^IBRITE^IKEGG_TC^ICAZy^IBiGG_Reaction".sp
    ...: lit()
egg = {}
     ...: for l in open("./emapper.tsv"):
     ...:          if l.startswith("#"):
     ...:                    continue
     ...:
     ...:          r = l.split("\t")
     ...:          x = {k:r[i] for i,k in enumerate(cols)}
     ...:          x["desc"] = " ".join( r[21 :]) .strip()
     ...:          egg[r[0]] = x
    ...: #prots = bpio.to_dict(bpio.parse("./tricho.faa","fasta"))
    ...: hits_ref = defaultdict(lambda:[])
    ...: for x in open("tricho_ref.tbl"):
    ...:     q,h = x.split()[:2]
    ...:     hits_ref[q].append(h.split("|")[2])
    ...: hits_ref = dict(hits_ref)
    ...: with  open("tricho.ann.gff3","w") as h:#, open("tricho.pw.gb","w") as h2:
    ...:     records = []
    ...:     for contig in tqdm(GFF.parse("./tricho.gff")):
    ...:         for f in contig.features:
    ...:             f.qualifiers["locus_tag"] = f.id
    ...:             for sf in f.sub_features:
    ...:                 qs = sf.qualifiers
    ...:                 qs["locus_tag"] = [sf.id]
    ...:                 if sf.id in hits_ref:
    ...:                     qs["hits_ref"] = hits_ref[sf.id]
    ...:                 if sf.id in prots:
    ...:                     qs["translation"] = [str(prots[sf.id].seq)]
    ...:                 if sf.id in egg:
    ...:                     qs["product"] = [egg[sf.id]["Preferred_name"]]
    ...:                     assert "EC_number" not in qs
    ...:                     if egg[sf.id]["EC"].strip():
    ...:
    ...:                         qs["EC_number"]  = egg[sf.id]["EC"].split(",")
    ...:
    ...:
    ...:                     if  egg[sf.id]["GOs"].strip():
    ...:                         qs["GO"] = egg[sf.id]["GOs"].split(",")
    ...:                     qs["seed_eggNOG_ortholog"] = egg[sf.id]["seed_eggNOG_ortholog"]
    ...:                     for x in ["KEGG_ko","KEGG_TC","KEGG_Reaction","KEGG_Pathway","KEGG_Module","CAZy","BiGG_Reaction","BRITE"]:
    ...:                         if egg[sf.id][x]:
    ...:                             qs[x] = egg[sf.id][x].split(",")
    ...:
    ...:                 if sf.id in interprot:
    ...:                     qs["note"] =  " | ".join( list(set([x.dbdesc for x in interprot[sf.id] if x.dbdesc])) )
    ...:                 if sf.id in antismash:
    ...:                     remove = "codon_start  locus_tag  product   protein_id  transl_table  translation".split()
    ...:                     for q in antismash[sf.id]:
    ...:                         if q not in remove:
    ...:                             qs["antismash_" + q ] = list(set(antismash[sf.id][q]))
    ...:         records.append(contig)
    ...:     GFF.write(records,h)

"""

"""
import pandas as pd                                                                                                                                                                                                                  
    ...: #todo = bpio.to_dict( GFF.parse("./tricho.gff"))                                                                                                                                                                                     
    ...: #interprot2 = pd.read_csv("interprot.tsv",sep="\t",names=["q","x1","length","db","dbid","dbdesc","start","end", "x2","x3","date"])                                                                                                   
    ...: #interprot2 = interprot2.fillna("")                                                                                                                                                                                                  
    ...: #interprot = defaultdict(lambda:[])                                                                                                                                                                                                  
    ...: #for _,x in interprot2.iterrows():                       
    ...: #    interprot[x.q].append(x)                            
    ...: #egg = {}                                   
    ...: #antismash = {}                                          
    ...: #with open("")                                           
    ...: #for x in bpio.parse("antismash.gbk","gb"):              
    ...: #    for f in x.features:                                
    ...: #        if "protein_id" in f.qualifiers:                                                  
    ...: #            antismash[f.qualifiers["protein_id"][0].replace("ncbi_","")] = f.qualifiers   
    ...:#for l in open("./emapper.tsv"):                         
    ...: #    if l.startswith("#"):                               
    ...: #        continue                                                                          
    ...: #    r = l.split("\t")                                                                     
    ...: #    x = {k:r[i] for i,k in enumerate(cols)}                                               
    ...: #    egg[r[0]] = x                                                                         
    ...: #prots = bpio.to_dict(bpio.parse("./tricho.faa","fasta"))                                  
    ...: hits_ref = defaultdict(lambda:[])                                                          
    ...: for x in open("tricho_ref.tbl"):                                                           
    ...:     q,h = x.split()[:2]                                                                    
    ...:     hits_ref[q].append(h.split("|")[2])                                                    
    ...: hits_ref = dict(hits_ref)                         
    ...: from Bio.Alphabet import generic_dna                                                       
    ...: genoma = bpio.to_dict(bpio.parse("./Trichoderma_atrovoridae_unknown.scaffolds.fa","fasta"))
    ...: with   open("tricho.pw.gb","w") as h2:      
    ...:     #records = []                                 
    ...:     for contig in tqdm(GFF.parse("./tricho.gff")):
    ...:         contig.seq = genoma[contig.id].seq  
    ...:         contig.seq.alphabet = generic_dna   
    ...:         features = []                               
    ...:         for f in contig.features:                   
    ...:                                                                
    ...:             f.qualifiers["locus_tag"] = f.id                   
    ...:             for sf in f.sub_features:                             
    ...:                 features.append(sf)                               
    ...:                 if sf.type == "mRNA":                                  
    ...:                     sf.type = "CDS"    
        ...:                     cds = [x.location for x in sf.subfeatures if x.type == "CDS"]
    ...:                     if len(cds) > 1:
    ...:                         sf.location = CompoundLocation( cds )
                                
    ...:                 qs = sf.qualifiers                                     
    ...:                 qs["locus_tag"] = [sf.id]                              
    ...:                 if sf.id in hits_ref:                                     
    ...:                     qs["hits_ref"] = hits_ref[sf.id]                      
    ...:                 if sf.id in prots:                                                                                        
    ...:                     qs["translation"] = [str(prots[sf.id].seq)]                                                           
    ...:                 if sf.id in egg:                                                                                          
    ...:                     qs["product"] = [egg[sf.id]["Preferred_name"]]                                                        
    ...:                                                                                                                           
    ...:                     if egg[sf.id]["EC"].strip():                                                                          
    ...:                         qs["EC_number"]  = egg[sf.id]["EC"].split(",")                                                         
    ...:                                                                                                                                
    ...:                                                                                                                                
    ...:                     if  egg[sf.id]["GOs"].strip():                                                                             
    ...:                                                                                                                                
    ...:                         qs["go_component"] = [process_go_cc(x) for x in  egg[sf.id]["GOs"].split(",") if process_go_cc(x)]     
    ...:                         qs["go_function"] = [process_go_mf(x) for x in  egg[sf.id]["GOs"].split(",") if process_go_mf(x)]      
    ...:                         qs["go_process"] = [process_go_bp(x) for x in  egg[sf.id]["GOs"].split(",") if process_go_bp(x)]       
    ...:                                                                                                                                
    ...:                                                                                                                
    ...:                     #qs["seed_eggNOG_ortholog"] = egg[sf.id]["seed_eggNOG_ortholog"]                                           
    ...:                     for x in ["KEGG_ko","KEGG_TC","KEGG_Reaction","KEGG_Pathway","KEGG_Module","CAZy","BiGG_Reaction","BRITE"]:
    ...:                         if egg[sf.id][x]:                                                                      
    ...:                             qs[x] = egg[sf.id][x].split(",")                                                   
    ...:                                                                                                                
    ...:                 if sf.id in interprot:                                                                         
    ...:                     qs["note"] =  " | ".join( list(set([x.dbdesc for x in interprot[sf.id] if x.dbdesc])) )
    ...:                 #if sf.id in antismash:                                                                        
    ...:                 #    remove = "codon_start  locus_tag  product   protein_id  transl_table  translation".split()
    ...:                 #    for q in antismash[sf.id]:
    ...:                 #        if q not in remove:                                       
    ...:                 #            qs["antismash_" + q ] = list(set(antismash[sf.id][q]))
    ...:         #records.append(contig)   
    ...:         contig.features = features
    ...:         bpio.write(contig,h2,"gb")
    ...: #    GFF.write(records,h)

"""