from setuptools import setup,find_packages

__version__ = "Undefined"
for line in open('SNDG/__init__.py'):
    if (line.startswith('__version__')):
        exec(line.strip())

setup(
    name='sndg-bio',
    packages=find_packages(),
    version=__version__,
    scripts=["scripts/structurome.py", "scripts/vcf2aln.py", "scripts/render_tree.py"],
    description="My bioinformatic scripts",
    author='Ezequiel Sosa - SNDG',
    install_requires=["tqdm", "bcbio-gff", "biopython", "goatools","networkx","python-libsbml"],
    author_email='ezequieljsosa@gmail.com',
    url='https://github.com/ezequieljsosa/sndg-bio',
    download_url='https://github.com/ezequieljsosa/sndg-bio/archive/0.1.tar.gz',
    keywords=['bioinformatics', 'sequence', 'example'],
    classifiers=['Programming Language :: Python', 'Topic :: Scientific/Engineering :: Bio-Informatics',
                 'Intended Audience :: Science/Research', ],
)

# rm -r dist/;python setup.py sdist bdist_wheel;twine upload dist/*
# ~/.pypirc
# [pypi]
# repository:https://upload.pypi.org/legacy/
# username = xxxxx
# password = xxxxx
