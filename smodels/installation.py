#!/usr/bin/env python

"""
.. module:: installation
   :synopsis: a module for returning installation paths and version numbers.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from __future__ import print_function
import sys
import os


def installDirectory():
    """
    Return the software installation directory, by looking at location of this
    method.
    
    """
    #path = os.path.abspath(os.path.join(os.path.realpath(__file__), os.pardir))
    path = os.path.abspath(__file__)
    path = os.path.abspath(os.path.join(path, '../..'))
    #path = path.replace("EGG-INFO/scripts/smodels-config", "")
    #path = path.replace("installation.py", "")
    return path + "/"


def pythonDirectory():
    """
    Return the python installation directory, by looking at location of this
    method. Same as installDirectory(), but trailing "smodels/" removed.
    
    """
    path = installDirectory()
    # path = path.replace("/smodels/", "/")
    return path


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


def banner():
    """
    Print SModelS banner.
    
    """
    f = open(installDirectory() + "BANNER")
    lines = f.readlines()
    f.close()
    return "".join(lines)


def printHelp():
    """
    Print usage information of this module.
    
    """
    print(sys.argv[0] + " [--help] [--installdir] [--pythondir]:")
    print("--help: show this message")
    print("--installdir: print SModelS installation directory")
    print("--pythondir: print SModelS python path")
    sys.exit(0)


if __name__ == "__main__":
    # print( banner() )
    if len(sys.argv) < 2:
        printHelp()
    for i in sys.argv[1:]:
        if i == "--help":
            printHelp()
        if i == "--installdir":
            print(installDirectory())
            sys.exit(0)
        if i == "--pythondir":
            print(pythonDirectory())
            sys.exit(0)
    printHelp()
