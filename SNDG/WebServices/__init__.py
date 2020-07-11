"""
Utilities for file download
"""

import hashlib
import logging
import os


from SNDG import execute

_log = logging.getLogger(__name__)

PROXIES = {}
def proxy_vars():
    vars = ""
    for k,v in PROXIES.items():
        vars= k + "=" + v
    return vars
def OvewriteFileException(Exception):
    pass


def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def md5_equal(fname, known_md5):
    match = md5(fname) == known_md5
    if not match:
        _log.error("error in %s checksum" % fname)
    return match


def download_file(complete_url, target, ovewrite=False, retries=3,timeout=20):
    if not target.strip():
        target = "./"
    if not os.path.exists(os.path.dirname(os.path.abspath(target))):
        raise FileNotFoundError("%s does not exists" % os.path.dirname(target))
    if os.path.exists(target) and not ovewrite:
        raise OvewriteFileException("%s already exists" % target)

    execute(f'wget  --timeout={timeout} --tries={retries} -O "{target}" "{complete_url}"')
