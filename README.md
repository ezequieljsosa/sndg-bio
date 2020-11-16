# sndg-bio
Bioinformatic tools

## Requirements
* Linux / tested on Ubuntu distribution 
* Python > 3.7 (sudo apt install python3.7)
* Python Pip (sudo apt install python3-pip)
* Docker (https://docs.docker.com/engine/install/ubuntu/)
* Virtual-env (pip3 install virtualenv) -> Recommended 

## Instalation
### From source

```console
git clone https://github.com/ezequieljsosa/sndg-bio.git
cd sndg-bio
virtual-env -p python3.7 env
echo export PYTHONPATH=$PYTHONPATH:$(pwd) >> env/bin/activate
source env/bin/activate
pip install -f requirements.txt
# OPTIONAL: set python env as default
echo source $(pwd)/env/bin/activate >> ${HOME}/.bashrc
```

## Commands

### Best bidirectional hits
Runs BBH algorithm using diamond-aligner
```console
# Help python3 -m "SNDG.Sequence.BBH" -h
python3 -m "SNDG.Sequence.BBH" -f1 proteome1.faa -f2 proteome2.faa > bbhs.txt
```





