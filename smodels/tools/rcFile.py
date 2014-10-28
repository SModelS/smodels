#!/usr/bin/env python

"""
.. module:: rcFile
   :synopsis: When imported, ~/.smodelsrc is parsed.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import os
import logging

logger = logging.getLogger(__name__)


def yesno(B):
    """
    TODO: write docstring
    
    """
    if B:
        return "yes"
    return "no"


def parseRCFile():
    """
    TODO: write docstring
    
    """
    rcfile = os.path.expanduser("~") + "/.smodelsrc"
    exists = os.path.exists(rcfile)
    if exists:
        execfile(rcfile)
        return True
    return False


if __name__ == "__main__":
    """
    Check if there is a smodelsrc file.
    
    """
    T = parseRCFile()
    logger.info("Checking if a ~/.smodelsrc file exists:")
    if not T:
        logger.warning("no ~/.smodelsrc file found.")
    else:
        logger.info("found ~/.smodelsrc file.")
