from setuptools import setup

with open("README.md") as handle:
    readme_md = handle.read().replace("\n","")

__version__ = "Undefined"
for line in open('SNDG/__init__.py'):
    if (line.startswith('__version__')):
        exec(line.strip())

setup(
    name = 'sndg-bio',
    packages = ['SNDG','SNDG.Structure','SNDG.Comparative','SNDG.WebServices'],
    version = __version__,
    scripts = ["scripts/structurome.py","scripts/vcf2aln.py"],
    description = readme_md,
    author = 'Ezequiel Sosa - SNDG',
    install_requires=["tqdm","bcbio-gff","biopython","goatools"],
    author_email = 'ezequieljsosa@gmail.com',
    url = 'https://github.com/ezequieljsosa/sndg-bio',
    download_url = 'https://github.com/ezequieljsosa/sndg-bio/archive/0.1.tar.gz',
    keywords = ['bioinformatics', 'sequence', 'example'],
    classifiers = [ 'Programming Language :: Python','Topic :: Scientific/Engineering :: Bio-Informatics','Intended Audience :: Science/Research',],
)

# python setup.py sdist bdist_wheel
# twine upload dist/*
# ~/.pypirc
# [pypi]
# repository:https://upload.pypi.org/legacy/
# username = xxxxx
# password = xxxxx
