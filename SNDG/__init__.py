import logging
import os
import subprocess as sp
import warnings
from subprocess import CalledProcessError
import sys
import shutil
import fileinput

from Bio import BiopythonWarning, BiopythonExperimentalWarning, BiopythonParserWarning

warnings.simplefilter('ignore', BiopythonWarning)
warnings.simplefilter('ignore', BiopythonExperimentalWarning)
warnings.simplefilter('ignore', BiopythonParserWarning)

log_format = "%(asctime)s - %(name)s - %(lineno)d - %(levelname)s - %(message)s"

__version__ = '0.1.57'

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
    "bgzip": "biocontainers/htslib:v1.2.1_cv3",
    "bedtools": "biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1",
    "bcftools": "biocontainers/bcftools:v1.9-1-deb_cv1",
    "spades": "staphb/spades:latest",
    "spades.py": "staphb/spades:latest",
    "bwa": "biocontainers/bwa:v0.7.17_cv1",
    "blast": "ncbi/blast:latest",
    "makeblastdb": "ncbi/blast:latest",
    "mafft": "biocontainers/mafft:v7.407-2-deb_cv1",
    "diamond-aligner": "biocontainers/diamond-aligner:v0.9.24dfsg-1-deb_cv1",
    "gatk": "broadinstitute/gatk:4.1.8.0"
}

import shlex


def _get_img(cmd):
    parts = shlex.split(cmd)
    img = [v for k, v in DOCKER_MAPPINGS.items() if k.startswith(parts[0])]
    return img[0] if img else None


def docker_wrap_command(cmd):
    if ";" in cmd:
        _log.warning("command containing ';', docker wrapping may fail")
    if "|" in cmd:
        for cmd2 in cmd.split("|"):
            img = _get_img(cmd2)
            if img:
                break
    else:
        img = _get_img(cmd)
    if not img:
        raise NotImplementedError(f"docker mapping for command {cmd} not found")

    parts = shlex.split(cmd, posix=False)
    new_parts = []
    mappings = []
    for p in parts:
        p = p[1:] if p[0] in ["'", '"'] else p
        p = p[:-1] if p[-1] in ["'", '"'] else p
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
    new_part = new_part.replace('"', '\\"')
    return f'docker run -u $(id -u):$(id -g) --rm -v $PWD:/out -w /out {mappings_str}  {img} bash -c "{new_part}"'


DEFAULT_SNDG_EXEC_MODE = "path_docker"


def execute(cmd, retcodes=[0], stdout=sys.stdout, stderr=sys.stderr,
            exec_mode=None):
    cmd = cmd.strip()
    exec_mode = exec_mode if exec_mode else os.environ.get("SNDG_EXEC_MODE", DEFAULT_SNDG_EXEC_MODE)
    if (exec_mode == "print"):
        print(cmd)
        return 0
    elif (exec_mode == "path_docker") and not cmd.startswith("docker"):
        parts = shlex.split(cmd, posix=False)
        if not shutil.which(parts[0]):
            _log.debug(f"command not detected, wrapping in docker:{cmd}")
            cmd = docker_wrap_command(cmd)

    elif (exec_mode == "docker") and not cmd.startswith("docker"):
        _log.debug(f"wrapping in docker:{cmd}")
        cmd = docker_wrap_command(cmd)
    try:
        _log.debug(cmd)
        process = sp.Popen(cmd, shell=True, stdout=stdout, stderr=stderr)
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

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return str(self.__dict__)


def grouper(iterable, n, fillvalue=None):
    from itertools import izip_longest
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return izip_longest(*args, fillvalue=fillvalue)


def arg_file_iter(arg_input):
    if arg_input == "-":
        file_iter = fileinput.input(arg_input)
    elif len(arg_input) > 1:
        file_iter = arg_input
    else:
        with open(arg_input, "r") as file:
            first_line = file.readline()
        if os.path.exists(first_line):
            with open(arg_input, "r") as file:
                file_iter = [x.strip() for x in file.readlines()]
        else:
            file_iter = [arg_input]
    return file_iter

