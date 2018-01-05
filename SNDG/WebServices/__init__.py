"""
Utilities for file download
"""

import os
import logging

from SNDG import execute

_log = logging.getLogger(__name__)


def OvewriteFileException(Exception):
    pass


def download_file(complete_url,target,ovewrite=False,retries=3):

    if not os.path.exists( os.path.dirname(target)):
        raise Exception("%s does not exists" % os.path.dirname(target))
    if os.path.exists(target) and not ovewrite:
        raise OvewriteFileException("%s already exists" % target)

    execute("wget -q --tries={retries} -O {target} {url}",
            url=complete_url,retries=retries, target=target)