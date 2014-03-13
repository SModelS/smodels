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
