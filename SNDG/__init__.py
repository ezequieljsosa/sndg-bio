import logging
import subprocess as sp

log_format = "%(asctime)s - %(name)s - %(lineno)d - %(levelname)s - %(message)s"

__version__ = '0.1'

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
    sp.call(cmd, shell=True)


class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)
