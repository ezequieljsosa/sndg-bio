In [103]: import subprocess as sp
...: sys.path.append("/home/eze/workspace/sndg-bio/")
...: from SNDG import mkdir
...: results = {}
...: for assembly in glob("groups/arg_chi/assemblies/*.fasta"):
    ...:     name = assembly.split("/")[-1].split(".fasta")[0]
...:     mkdir(f"strains/{name}/blast_mec/")
...:     assert os.path.exists(f"strains/{name}/blast_mec/")
...:     sp.check_output(f"blastn -qcov_hsp_perc 10 -perc_identity 90 -outfmt '6 qseqid sseqid pident length qcovs pident qstart qend sstart send' -evalue 1e-6 -out strains/{name}/blast_mec/result.t
...: bl -db  {assembly} -query ./curated/mec_sequences/mecs.fna",shell=True)
    ...:     results[name] = pd.read_table(f"./strains/{name}/blast_mec/result.tbl",names='qseqid sseqid pident length qcovs qstart qend sstart send'.split())
    ...:
...:
...:
...:


In [103]: from itertools import product
...: def aligned_subcontigs(df):
    ...:     finished = False
...:     hits = []
...:     print(df)
...:     df_tmp = df.copy()

...:     df_tmp["start"] = [r.sstart if r.sstart < r.send else r.send for _,r in df_tmp.iterrows()]
...:     df_tmp["end"] = [r.send if r.sstart < r.send else r.sstart for _,r in df_tmp.iterrows()]
...:     df_final = df
...:     while not finished:
    ...:         best_contig = df_tmp.groupby("qseqid").agg({"length":"sum"}).sort_values("length").iloc[-1].name
...:         df2 = [x.to_dict() for _,x in df_tmp[df_tmp.qseqid == best_contig].iterrows()]
...:         inicial = len(df2)
...:         done = False
...:         print(pd.DataFrame(df2))
...:         while not done:
    ...:             done = True
...:             for x,y in [(x,y) for x,y in product(df2,df2) if x != y]:
    ...:                 range1 = set([ x["sseqid"] + "|" + str(j) for j in range(x["start"],x["end"])])
...:                 range2 = set([ y["sseqid"] + "|" + str(j) for j in range(y["start"],y["end"])])
...:                 if len(range1 & range2) > 0:
    ...:                     if len(range1 & range2) > 500:
    ...:                         print(df2)
...:                         raise Exception("afklsddf")
...:                     df2 = [z for z in df2 if z not in [x,y] ]
...:                     r = x
...:                     r["start"] = min([x["start"],y["start"]])
...:                     r["end"] = max([x["end"],y["end"]])
...:                     df2 = df2.append(r)
...:                     df2 = pd.DataFrame(df2)
...:                     done = False
...:                     break
...:
...:         df2 = pd.DataFrame(df2)
...:         print(df2)
...:         break
...:         assert len(df2) <= inicial, [len(df2), inicial]
...:         df_final = []
...:         df3 = df_tmp[df_tmp.qseqid != best_contig]
...:         for _,r1 in df3.iterrows():
    ...:            agregar  = True
...:            for _,r2 in df2.iterrows():
    ...:                range1 = set([ r1["sseqid"] + "|" + str(j) for j in range(r1["start"],r1["end"])])
...:                range2 = set([ r2["sseqid"] + "|" + str(j) for j in range(r2["start"],r2["end"])])
...:                if len(range1 & range2) > 500:
    ...:                    agregar = False
...:            if agregar:
    ...:                    df_final.append(r1)
...:                    break
...:         assert len(df_final) <= len(df3) , [len(df_final) , len(df3)]
...:         xx = pd.DataFrame(df_final)
...:         df_final = pd.concat( [xx , df3] )
...:         print(len(df_final),len(df_tmp))
...:
...:         assert len(df_tmp) >= len(df_final), [len(df_final) , len(df_tmp)]
...:         if len(df_final) == len(df_tmp):
    ...:             finished = True
...:         else:
...:             df_tmp = df_final
...:     return df_final
...: df = aligned_subcontigs(results["0037"])