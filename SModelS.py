#!/usr/bin/env python

"""
.. module:: SModelS
   :synopsis: Intended as a potential main entry point, currently just for
   returning the SModelS version number.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from __future__ import print_function
import sys
import os
import inspect

def installDirectory():
    """
    Return the software installation directory, by looking at location of this
    method.
    
    """
    ret = os.path.realpath(inspect.getabsfile(installDirectory))
    ret = ret.replace("SModelS.py", "")
    return ret


def version(astuple=False):
    """
    Print version number of the SModelS framework.
    
    """
    f = open("%s/version" % installDirectory())
    l = f.readline()
    f.close()
    l = l.replace("\n", "")
    l.strip()
    if not astuple:
        return l
    a, b = l.split(".")
    return (int(a), int(b))
        

def license():
    """
    Print license information of the SModelS framework.
    
    """
    f = open(installDirectory() + "COPYING")
    lines = f.readlines()
    f.close()
    return "".join(lines)


def printHelp():
    """
    Print usage information of this module.
    
    """
    print(sys.argv[0] + ": --help --installdir")
    print("--help: show this message")
    print("--installdir: print SModelS installation directory")
    sys.exit(0)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        printHelp()
    for i in sys.argv[1:]:
        if i == "--help":
            printHelp()
        if i == "--installdir": 
            print(installDirectory())
            sys.exit(0)
    printHelp()
