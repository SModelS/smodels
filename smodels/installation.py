#!/usr/bin/env python

"""
.. module:: installation
   :synopsis: a module for returning installation paths and version numbers.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from __future__ import print_function
from smodels.tools.smodelsLogging import logger
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


def authors():
    """ return the author list, taken from BANNER """
    copying_file = open('%s/smodels/share/BANNER' % installDirectory(), 'r')
    lines = copying_file.readlines()
    copying_file.close()
    authors = ""
    start_parsing=False
    for line in lines:
        if "Copyright" in line:
            start_parsing = True
        if not start_parsing: continue
        to_add = line.replace ( " <smodels-users@lists.oeaw.ac.at>","" )
        to_add = to_add.replace ( "Copyright (C) ","").replace ( "\n", "" )
        if to_add[:5]=="2012-":
            to_add = to_add[10:]
        if len(authors)>0:
            authors+=" "
        authors += to_add
    return authors

def _toTuple_ ( ver ):
    """ convert version string to tuple """
    a = ver.replace(" ",".",1).split(".")
    for ctr,el in enumerate(a):
        try:
            a[ctr]=int(el)
        except ValueError:
            a[ctr]=el
    b=[]
    for i in a:
        found=False
        for pf in [ "rc", "post", "pre" ]:
            if type(i)==str and pf in i:
                found=True
                minor = i[:i.find(pf)]
                try:
                    minor = int(minor)
                except:
                    pass
                b.append ( minor )
                b.append ( i[i.find(pf):] )
                continue
        if not found:
            b.append ( i )
    return tuple(b)

def requirements():
    ret=[]
    f = open("%s/smodels/share/requirements.txt" % installDirectory())
    lines=f.readlines()
    for l in lines: ret.append ( l.strip() )
    f.close()
    return ret

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
    return _toTuple_ ( l )


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
    f = open(installDirectory() + "/smodels/share/BANNER")
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

def fixpermissions():
    """ make sure that all filepermissions are such that
        we can compile the wrappers for pythia and nllfast. """
    import os, glob
    Dir = "%ssmodels/lib/" % installDirectory()
    try:
        Dirs = [ "%spythia6" % Dir, "%spythia8" % Dir ]
        Dirs += glob.glob("%snllfast/nllfast-*" % Dir )
        Dirs += glob.glob("%spythia8/xml.doc" % Dir )
        for p in Dirs:
            logger.debug ( "chmod 777 %s" % (p) )
            os.chmod ( p, 0o777 )
    except Exception as e:
        print ( "chmod failed (permission error). Please try as root, i.e.:" )
        print ( "sudo smodelsTools.py fixpermissions" )

def main():
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
        if i in [ "--requirements", "-r" ]:
            print(requirements())
            sys.exit(0)
        if i in [ "--help", "-h" ]:
            printHelp()
            sys.exit(0)
        if i in [ "--license", "--copyright", "-c" ]:
            print(license())
            sys.exit(0)
    print("Error: cannot parse %s.\n" % i )
    printHelp()

if __name__ == "__main__":
    main()
