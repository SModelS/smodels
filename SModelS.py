#!/usr/bin/env python

"""
.. module:: SModelS
    :synopsis: Intended as a potential main entry point, currently just for
               returning the SModelS version number.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

def installDirectory():
    """ return the software installation directory, by looking at location of this method """
    import os, inspect
    ret=os.path.realpath ( inspect.getabsfile(installDirectory) )
    ret=ret.replace("SModelS.py","")
    return ret

def version( astuple=False, addCodeName=True ):
    """ prints out version number of SModelS framework """
    f=open("%s/version" % installDirectory() )
    l=f.readline()
    f.close()
    l=l.replace("\n","")
    l.strip()
    if not astuple:
        if addCodeName: return l
        p=l.find("/")
        if p>-1: return l[:p]
    T,C=l.split("/")
    A,B=T.split(".")
    if addCodeName:
        return (int(A),int(B),C.strip())
    else:
        return (int(A),int(B))
        

def license():
    f=open(installDirectory()+"COPYING")
    lines=f.readlines()
    f.close()
    return "".join(lines)

def help():
  import sys
  print sys.argv[0]+": --help --installdir"
  print "--help: show this message"
  print "--installdir: print SModelS installation directory"
  sys.exit(0)


if __name__ == "__main__":
    import sys
    if len(sys.argv)<2: help()
    for i in sys.argv[1:]:
        if i=="--help": help()
        if i=="--installdir": print installDirectory()

