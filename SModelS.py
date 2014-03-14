#!/usr/bin/env python

"""
.. module:: SModelS
    :synopsis: Intended as a potential main entry point, currently just for
               returning the SModelS version number.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

def installdir():
    """ return the software installation directory, by looking at location of this method """
    import os, inspect
    ret=os.path.realpath ( inspect.getabsfile(installdir) )
    ret=ret.replace("SModelS.py","")
    return ret

def version( astuple=False ):
    """ prints out version number of SModelS framework """
    f=open("%s/version" % installdir() )
    l=f.readline()
    f.close()
    l=l.replace("\n","")
    l.strip()
    if not astuple: return l
    T,C=l.split("/")
    A,B=T.split(".")
    return (int(A),int(B),C.strip())

def license():
    f=open(installdir()+"COPYING")
    lines=f.readlines()
    f.close()
    return "".join(lines)

def installdir():
  import os, inspect
  ret=os.path.realpath ( inspect.getabsfile(installdir) )
  ret=ret.replace("bin/smodels-config","").replace("SModelS.py","")
  print ret

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
        if i=="--installdir": installdir()

