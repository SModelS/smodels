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
    path = os.path.abspath(os.path.realpath(__file__))
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
    f = open("%s/smodels/version" % installDirectory())
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
    Returns SModelS banner.
    
    """
    f = open(installDirectory() + "BANNER")
    lines = f.readlines()
    f.close()
    return "".join(lines)

def printHelp():
    """
    Print usage information of this module.
    
    """
    print("Usage: " + sys.argv[0] + " [--help|-h] [--installdir|-i] [--pythondir|-p]")
    print("                      [--version|-v] [--banner|-b] [--license|--copyright|-c]:")
    print("--help:       show this message")
    print("--installdir: print SModelS installation directory")
    print("--pythondir:  print SModelS python path")
    print("--version:    print SModelS version number")
    print("--banner:     print SModelS banner")
    print("--copyright:  print SModelS copyright")
    sys.exit(0)


if __name__ == "__main__":
    # print( banner() )
    if len(sys.argv) < 2:
        printHelp()
    for i in sys.argv[1:]:
        if i in [ "--installdir", "-i" ]:
            print(installDirectory())
            sys.exit(0)
        if i in [ "--pythondir", "-p" ]:
            print(pythonDirectory())
            sys.exit(0)
        if i in [ "--version", "-v" ]:
            print(version())
            sys.exit(0)
        if i in [ "--banner", "-b" ]:
            print(banner())
            sys.exit(0)
        if i in [ "--help", "-h" ]:
            printHelp()
            sys.exit(0)
        if i in [ "--license", "--copyright", "-c" ]:
            print(license())
            sys.exit(0)
    print("Error: cannot parse %s.\n" % i )
    printHelp()
