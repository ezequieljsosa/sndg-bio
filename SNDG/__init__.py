import logging
import os
import subprocess as sp
import warnings
from subprocess import CalledProcessError
import sys
import time

from Bio import BiopythonWarning, BiopythonExperimentalWarning, BiopythonParserWarning

warnings.simplefilter('ignore', BiopythonWarning)
warnings.simplefilter('ignore', BiopythonExperimentalWarning)
warnings.simplefilter('ignore', BiopythonParserWarning)

log_format = "%(asctime)s - %(name)s - %(lineno)d - %(levelname)s - %(message)s"

__version__ = '0.1.27'

_log = logging.getLogger(__name__)

test_execution = []


def init_log(log_file_path=None, rootloglevel=logging.DEBUG):
    default_formatter = logging.Formatter(log_format)
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(default_formatter)
    root = logging.getLogger()

    if log_file_path:
        fh = logging.FileHandler(log_file_path)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(default_formatter)
        root.addHandler(fh)

    root.addHandler(console_handler)
    root.setLevel(rootloglevel)


DOCKER_MAPPINGS = {
    "samtools": "biocontainers/samtools:v1.9-4-deb_cv1",
    "bedtools": "biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1",
    "bcftools": "biocontainers/bcftools:v1.9-1-deb_cv1",
    "spades": "staphb/spades:latest",
    "spades.py": "staphb/spades:latest",
    "bwa": "biocontainers/bwa:v0.7.17_cv1",
    "blast": "ncbi/blast:latest",
    "makeblastdb": "ncbi/blast:latest"
}

import shlex


def docker_wrap_command(cmd):
    if ";" in cmd:
        _log.warning("command containing ';', docker wrapping may fail")
    parts = shlex.split(cmd)
    img = [v for k, v in DOCKER_MAPPINGS.items() if k.startswith(parts[0])]
    if not img:
        raise NotImplemented(f"docker mapping for command {parts[0]} not found")
    else:
        img = img[0]

    new_parts = []
    mappings = []
    for p in parts:
        np = p

        if p.startswith("./") or p.startswith("../") or p.startswith("/"):
            np = os.path.abspath(p)
            if os.path.isdir(np):
                mappings.append(np)
            else:
                mappings.append(os.path.abspath(np + "/../"))
        elif p.startswith("~/"):
            np = os.path.expanduser(p)
            if os.path.isdir(np):
                mappings.append(np)
            else:
                mappings.append(os.path.abspath(np + "/../"))
        for arg_part in ["=/", "=./", "=../"]:
            if arg_part in p:
                parts2 = p.split(arg_part)
                corr = os.path.abspath(arg_part[1:] + parts2[1])
                np = parts2[0] + "=" + corr
                if os.path.isdir(corr):
                    mappings.append(corr)
                else:
                    mappings.append(os.path.abspath(corr + "/../"))

        new_parts.append(np)

    final_mappings = mappings.copy()
    change = True
    while change:
        change = False

        final_mappings = {m: 1 for m in final_mappings}
        for mapping1 in final_mappings.copy():
            for mapping2 in final_mappings.copy():
                if mapping1 != mapping2:
                    if (mapping1 in mapping2) and (mapping1 in final_mappings):
                        del final_mappings[mapping1]
                        change = True
                    elif (mapping2 in mapping1) and (mapping2 in final_mappings):
                        del final_mappings[mapping2]
                        change = True
        final_mappings = list(final_mappings)

    mappings_str = " ".join([f"-v {x}:{x}" for x in set(final_mappings)])

    new_part = " ".join([x if " " not in x else f'"{x}"' for x in new_parts])
    return f'docker run -u $(id -u):$(id -g) --rm -v $PWD:/out -w /out {mappings_str}  {img} {new_part}'


def execute(cmd_unformated, wd="./", retcodes=[0], docker_mode=False, **kargs):
    cmd = cmd_unformated.format(**kargs)

    try:
        process = sp.Popen(cmd, shell=True)
        process.communicate()
        return process.returncode

    except CalledProcessError as ex:
        if ex.returncode not in retcodes:
            raise
    return ex.returncode


def execute_from(cmd_unformated, workdir, **kargs):
    cwd = os.getcwd()
    try:
        os.chdir(workdir)
        execute(cmd_unformated, **kargs)
    finally:
        os.chdir(cwd)


def mkdir(dirpath):
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)


class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)


def grouper(iterable, n, fillvalue=None):
    from itertools import izip_longest
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return izip_longest(*args, fillvalue=fillvalue)
