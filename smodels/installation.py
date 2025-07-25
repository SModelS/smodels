#!/usr/bin/env python3

"""
.. module:: installation
   :synopsis: a module for returning installation paths and version numbers.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from __future__ import print_function
import sys
import os
from typing import Union, Tuple

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

def test_requirements():
    """ checks if all requirements are installed.
    Returns true if that is the case. """
    import importlib.util
    for i in requirements():
        pos = i.find(">" )
        lib = i
        if pos > -1:
            lib = i[:pos]
        found = importlib.util.find_spec( lib )
        if found == None:
            return False
    return True

def resolve_dependencies( as_user = True, verbose = True ):
    """ method that is meant to resolve the SModelS dependencies,
    via pip install --user. Warns you if pip cannot be found.
    :params as_user: if False, try system-wide install.
    """
    ck = test_requirements()
    if ck == True: ## nothing to be done
        print ( "all dependencies are met." )
        return None
    import subprocess
    req = f"{installDirectory()}/smodels/share/requirements.txt"
    try:
        import pip
    except (ModuleNotFoundError,ImportError) as e:
        if verbose:
            print ( "warning: pip not found. cannot install requirements. Maybe try easy_install pip" )
        return
        # sys.exit()
    cmd = [ sys.executable, '-m', 'pip', 'install', '--upgrade', '-r', req ]
    if as_user:
        o = subprocess.getoutput ( " ".join ( [ sys.executable, '-m', 'pip', 'install', '--user', '--dry-run', 'pyslha' ] ) )
        #o = subprocess.getoutput ( f"{p} install --user --dry-run pyslha" )
        if not "Can not perform a '--user' install" in o:
            cmd.insert( 4, '--user' )

    out = subprocess.call ( cmd )
    if out == 0:
        print ( "dependencies have been installed successfully." )
        return None
    else:
        print ( "an error has occurred when resolving the dependencies." )
        return -1

def cacheDirectory ( create=False, reportIfDefault=False ):
    """
    Returns the user's smodels cache directory, i.e. ~/.cache/smodels.
    :params create: if True, create the directory if it doesnt exist.
    :params reportIfDefault: if True, then report also if the default location
                             has been used
    :returns: cache dir. optionally returns also if the default cache dir has been
              used.
    """
    if "SMODELS_CACHEDIR" in os.environ:
        cacheDir = os.environ["SMODELS_CACHEDIR"]
        if create and not os.path.exists ( cacheDir ):
            os.mkdir ( cacheDir )
        if reportIfDefault:
            return cacheDir,False
        return cacheDir
    home = os.environ["HOME"]
    cacheDir = os.path.join ( home, ".cache" )
    if create and not os.path.exists ( cacheDir ):
        os.mkdir ( cacheDir )
    smodelsDir = os.path.join ( cacheDir, "smodels" )
    if create and not os.path.exists ( smodelsDir ):
        os.mkdir ( smodelsDir )
    if reportIfDefault:
        return smodelsDir,True
    return smodelsDir

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
    copying_file = open(f'{installDirectory()}/smodels/share/BANNER', 'r')
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

def requirements():
    ret=[]
    f = open( f"{installDirectory()}/smodels/share/requirements.txt")
    lines=f.readlines()
    for l in lines: ret.append ( l.strip() )
    f.close()
    return ret

def version( return_tuple : bool = False ) -> Union[str,Tuple]:
    """
    Returns version number of the SModelS framework.

    :param return_tuple: Return a tuple of (major,minor,patch,...) instead of string
    :returns: version, either as string or tuple
    """
    f = open( f"{installDirectory()}/smodels/version" )
    l = f.readline()
    f.close()
    l = l.replace("\n", "")
    l.strip()
    if not return_tuple:
        return l
    import re
    ret = re.split("\\.|-",l)
    for i,r in enumerate(ret):
        try:
            ret[i]=int(r)
        except ValueError as e:
            pass
    return tuple(ret)

def license():
    """
    Print license information of the SModelS framework.

    """
    f = open(installDirectory() + "smodels/COPYING")
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

def fixpermissions():
    """ make sure that all filepermissions are such that
        we can compile the wrappers for pythia and nllfast. """
    from smodels.base.smodelsLogging import logger
    import glob
    Dir = f"{installDirectory()}smodels/lib/"
    try:
        Dirs = [ f"{Dir}pythia6", f"{Dir}pythia8" ]
        Dirs += glob.glob(f"{Dir}nllfast/nllfast-*" )
        Dirs += glob.glob(f"{Dir}pythia8/xml.doc" )
        for p in Dirs:
            logger.debug ( f"chmod 777 {p}" )
            os.chmod ( p, 0o777 )
    except OSError:
        print ( "chmod failed (permission error). Please try as root, i.e.:" )
        print ( "sudo smodelsTools.py fixpermissions" )

__dbServer__ = "https://smodels.github.io/database"
__dblabels__ = [ "official", "latest", "fastlim", "backup", "superseded", "unittest",
                 "debug", "nonaggregated", "full_llhds", "unittestextra", None ]

def databasePath ( label ):
    """ construct the path to the database json file
    :param label: one of: official, latest, fastlim, backup
    :returns: URL, e.g. https://smodels.github.io/database/official
    """
    if not label in __dblabels__:
        from smodels.base.smodelsLogging import logger
        logger.warning ( f"cannot identify label {label}" )
        return label
    if label == None:
        label = "official"
    v=version().replace(".","")
    r=f"{__dbServer__}/{label}{v}"
    return r

def main():
    import argparse
    ap = argparse.ArgumentParser( description= "installation helper" )
    ap.add_argument( "-i", "--installdir", help="print SModelS installation directory", action="store_true" )
    ap.add_argument( "-p", "--pythondir", help="print SModelS python path",
                     action="store_true" )
    ap.add_argument( "-v", "--version", help="print SModelS version number",
                     action="store_true" )
    ap.add_argument( "-b", "--banner", help="print SModelS banner",
                     action="store_true" )
    ap.add_argument( "-R", "--resolve_dependencies", help="try to resolve the SModelS dependencies via pip",
                     action="store_true" )
    ap.add_argument( "-r", "--requirements",help="print SModelS python requirements",
                     action="store_true" )
    ap.add_argument( "-d", "--database",
                     help="print SModelS official database url for this release", action="store_true")
    ap.add_argument( "-t", "--test-database", help="print SModelS official unittest database url for this release", action="store_true" )
    ap.add_argument( "-c", "--copyright", "--license",
                     help="print SModelS copyright", action="store_true" )
    args = ap.parse_args()
    funcs = { "installdir": installDirectory, "pythondir": pythonDirectory,
              "version": version, "banner": banner, "requirements": requirements,
              "copyright": license, "resolve_dependencies": resolve_dependencies }
    dbpaths = { "database": "official", "test_database": "unittest" }
    for f,v in args.__dict__.items():
        if v:
            if "database" in f:
                t = databasePath ( dbpaths[f] )
                if t: print ( t )
            else:
                r = funcs[f]()
                if r != None:
                    print ( r )

if __name__ == "__main__":
    main()
