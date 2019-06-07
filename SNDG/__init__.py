import logging
import os
import subprocess as sp
import warnings
from subprocess import CalledProcessError

from Bio import BiopythonWarning, BiopythonExperimentalWarning, BiopythonParserWarning

warnings.simplefilter('ignore', BiopythonWarning)
warnings.simplefilter('ignore', BiopythonExperimentalWarning)
warnings.simplefilter('ignore', BiopythonParserWarning)

log_format = "%(asctime)s - %(name)s - %(lineno)d - %(levelname)s - %(message)s"

__version__ = '0.1.25'

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


def execute(cmd_unformated,wd="./",retcodes=[0], **kargs):
    cmd = "cd " + wd + ";" + cmd_unformated.format(**kargs)
    _log.debug( cmd)
    try:
        if not test_execution:
            sp.check_output(cmd, shell=True, stderr=sp.STDOUT)
        else:
            print(cmd)
        _log.debug(cmd + " -> OK")
    except CalledProcessError as ex:
        _log.warning(ex.output)
        if ex.returncode not in retcodes:
            raise


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
