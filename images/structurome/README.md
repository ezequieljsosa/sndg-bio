Structurome image has all the scripts needed to create model structures 
from a .fasta proteome.

It uses blast and psi blast on Uniref for the profile creation and alingment
and Modeller (https://salilab.org/modeller/) for homology modelling.

The script runs standalone, so it is necesary to download several databases
for it to run, a free space of 150 GB is required. The script createdb.py creates the directory structure, downloads 
the files and indexes and processes them. 
The required layout is the following:

base (any directory of your choosing, example: /home/sndg/data)   
├── uniprot  
│     ├── uniref90.fasta  
│     └── idmapping  
└── pdb  
    ├── entries.idx  
    ├── divided  
    │    ├── aa  
    │    ├── ax  
    │    │   ├── pdb1ax2.ent  
    │    │   ├── pdb7ax3.ent  
    │    │   └── ..  
    │    └── ..  
    └── processed  
         └── seqs_from_pdb.fasta      

To create the models, the image has several scripts. We recomend using 
structurome.py, which will process a fasta file to create models for 
complete sequences and domains 
(if no model for the complete sequence is avaliable).
Depending on the size of the proteome and the avaliable cpus, 
the process may take from days to weeks.
The "quick" option is avaliable, which will only produce models 
for proteins with high identity with a pdb chain. 


### TL;DR  
```{r, engine='bash', count_lines}
docker pull ezequieljsosa/structurome 
# with 150Gb avaliable, only the first time, it may take a few hours  
docker  run   --mount type=bind,source="/home/myuser/data",target="/data" ezequieljsosa/structurome createdb.py
docker  run   --mount type=bind,source="/home/myuser/data",target="/data",readonly --mount type=bind,source="/tmp/struct",target="/out" {valid_modeller_key} structurome.py /out/proteome.fasta  
```

# Development
Update distro:  
Change:  
```SNDG.__init__.__version__=x.y.z```
```{r, engine='bash', count_lines}
python setup.py bdist_wheel 
twine upload dist/*
```
