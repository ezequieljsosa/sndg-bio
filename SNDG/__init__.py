import os
import logging
import subprocess as sp
from subprocess import CalledProcessError

log_format = "%(asctime)s - %(name)s - %(lineno)d - %(levelname)s - %(message)s"

__version__ = '0.1.5'

_log = logging.getLogger(__name__)


def init_log(log_file_path=None):
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
    root.setLevel(logging.DEBUG)


def execute(cmd_unformated, **kargs):

    cmd = cmd_unformated.format(**kargs)
    _log.debug(cmd)
    try:
        sp.check_output(cmd, shell=True,stderr=sp.STDOUT)
        _log.debug(cmd  + " -> OK")
    except CalledProcessError as ex:
        _log.warn(ex.message)
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
