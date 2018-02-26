SNDG pannotator scripts are used to discover and annotate genes 
of a prokaryotic sequence.
The input is always a fasta, and it can be either an assembly or a proteome.
The generated outputs are a gff and a genebank file.

The script runs standalone, so it is necesary to download several databases
for it to run, a free space of 250 GB is required. The script createdb.py creates the directory structure, downloads 
the files and indexes and processes them. 
The required layout is the following:

base (any directory of your choosing, example: /home/sndg/data)  
├── uniprot  
│     ├── uniref  
│     │     ├── uniref90.fasta  
│     └── goa_uniprot_all.gpa.gz  
├── pfamtigrfam  
│     ├── pfamatigrfam.hmm  
│     ├── TIGRFAMS_GO_LINK  
│     └── INFO  
├── cog  
│     ├── myva  
│     └── whog  
└─── ec  
      └── PRIAM_MAR15   
            ├──  myva      
            └─── priam  

To create a container, 2 filesystem mappings are needed, 
one for the annotation data (must be mounted before running "createdb.py")
and other for the scripts output.  
To use the annotator one can use the annotate.sh script, 
which runs the entire pipeline or may choose
to run each script individually 

Ex:  
1) docker  run   --mount type=bind,source="/data",target="/data",readonly 
--mount type=bind,source="/tmp/docker",target="/out" 
pannotator annotate.sh /out/mygenome.fasta
In the machine, you must have the /tmp/docker/mygenome.fasta file, 
so it can be mapped to /out/mygenome.fasta in the container. 
This creates the files ncbi_mygenome.gff and mygenome.gb in /tmp/docker/

2)  docker  run   --mount type=bind,source="/data",target="/data",readonly 
--mount type=bind,source="/tmp/docker",target="/out"  pannotator bash  
Inside the container, you can run scripts /app/p_procariota/ from 1 to 8   

### TL;DR  
```{r, engine='bash', count_lines}
docker pull ezequieljsosa/pannotator 
# with 250Gb avaliable, only the first time, it may take a few hours  
docker  run   --mount type=bind,source="/home/myuser/data",target="/data" \
    ezequieljsosa/pannotator createdb.py
docker  run   --mount type=bind,source="/home/myuser/data",target="/data",readonly \
    --mount type=bind,source="/home/myuser/mygenome/",target="/out" \
    pannotator annotate.sh /out/mygenome.fasta
```

