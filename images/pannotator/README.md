SNDG pannotator scripts are used to discover and annotate genes 
of a prokaryotic sequence.
The input is always a fasta, and it can be either an assembly or a proteome.
The generated outputs are a gff and a genebank file.

The script runs standalone, so it is necesary to download several databases
for it to run, a free space of 150 GB is required. The script createdb.py creates the directory structure, downloads 
the files and indexes and processes them. 
The required layout is the following:

base (any directory of your choosing, example: /home/sndg/data)
├── uniprot
│     ├── uniref90.fasta
│     └── idmapping
├── pfamtigrfam
│     ├── pfamatigrfam.hmm
│     ├── TIGRFAMS_GO_LINK
│     └── INFO
├── cog
│     ├── myva
│     └─── whog
└─── priamrpsdb
      └─── priam
   
 


docker  run   --mount type=bind,source="/data",target="/data",readonly --mount type=bind,source="/tmp/docker",target="/out"  pannotator bash


