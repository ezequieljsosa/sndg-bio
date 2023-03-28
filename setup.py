from setuptools import setup, find_packages

__version__ = "2.0"

#REQUIREMENTS = [i.strip() for i in open("requirements.txt").readlines()]

setup(
    name='sndg-bio',
    packages=find_packages(),
    version=__version__,
    #scripts=["scripts/structurome.py", "scripts/vcf2aln.py", "scripts/render_tree.py"],
    description="My bioinformatic scripts",
    author='Ezequiel Sosa - SNDG',
    install_requires=["tqdm", "bcbio-gff", "biopython", "goatools", "dna_features_viewer"],
    author_email='ezequieljsosa@gmail.com',
    url='https://github.com/ezequieljsosa/sndg-bio',
    download_url='https://github.com/ezequieljsosa/sndg-bio/archive/0.1.tar.gz',
    keywords=['bioinformatics', 'sequence', 'example'],
    classifiers=['Programming Language :: Python', 'Topic :: Scientific/Engineering :: Bio-Informatics',
                 'Intended Audience :: Science/Research', ],
    extras_require={
        'mongo': ['pymongo', 'mongoengine'],
        'network': ['python-libsbml', "networkx"],
        'phylo': ['ete3'],
        'def': ["twine"]
    },
)

# rm -r dist/;python setup.py sdist bdist_wheel;twine upload dist/*
# ~/.pypirc
# [pypi]
# repository:https://upload.pypi.org/legacy/
# username = xxxxx
# password = xxxxx
