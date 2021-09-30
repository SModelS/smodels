#!/usr/bin/env python3

""" the installer script, fetches pythia8, explodes the tarball, compiles. """

import os, sys, shutil

def getVersion():
    """ obtain the pythia version we wish to use, stored in file 'pythiaversion' """
    if not os.path.exists ( "pythiaversion" ):
        print ( "[installer.py] error cannot determine pythiaversion." )
        sys.exit(-1)
    with open("pythiaversion","rt") as f:
        ver = f.read()
        ver = ver.strip()
        f.close()
    return ver

def rmTarball():
    """ remove cruft tarballs """
    ver = getVersion()
    tarball=f"pythia{ver}.tgz"
    if os.path.exists ( tarball ):
        os.unlink ( tarball )

def fetch():
    """ fetch tarball from pythia.org or smodels.github.io """
    ver = getVersion()
    tarball=f"pythia{ver}.tgz"
    if os.path.exists ( tarball ):
        size = os.stat ( tarball ).st_size
        if size > 15000000:
            print ( f"[installer.py] tarball {tarball} exists. Won't fetch." )
            return tarball
        else:
            rmTarball()
    import requests
    url=f"http://pythia.org/download/pythia{ver[:2]}/"
    print ( f"[installer.py] fetching {tarball} from {url}" )
    path = os.path.join ( url, tarball )
    # URL=http://home.thep.lu.se/~torbjorn/pythia8/$TARBALL
    #TARBALL="pythia${VER}_fixed.tgz"
    r = requests.get ( path, stream=True )
    if r.status_code != 200:
        print ( f"[installer.py] could not fetch tarball: {r.reason}." )
        rmTarball()
        url=f"http://smodels.github.io/downloads/tarballs/"
        print ( f"[installer.py] trying to fetch {tarball} from {url}" )
        path = os.path.join ( url, tarball )
        r = requests.get ( path, stream=True )
        if r.status_code != 200:
            print ( f"[installer.py] could not fetch tarball: {r.reason}." )
            rmTarball()
            sys.exit(-1)
    with open ( tarball, "wb" ) as f:
        shutil.copyfileobj( r.raw, f )
        f.close()

def unzip():
    """ explode the tarball """
    ver = getVersion()
    Dir = f"pythia{ver}"
    if os.path.exists ( Dir ):
        if os.path.isdir ( Dir ):
            print ( f"[installer.py] {Dir} exists, skip explosion of tarball." )
            return
        else:
            shutil.rmtree ( Dir )
    tarball=f"pythia{ver}.tgz"
    import tarfile
    print ( f"[installer.py] extracting to pythia{ver}" )
    with tarfile.open( tarball, 'r:gz') as f:
        f.extractall()

def getNCPUs():
    """ get the number of CPUs to compile pythia8 with """
    ncpus = 4
    try:
        from smodels.tools.runtime import nCPUs
        ncpus = nCPUS()
    except Exception as e:
        pass
    return ncpus

def checkPythia():
    """ check if a compiled pythia library already exists """
    ver = getVersion()
    afile = f"pythia{ver}/lib/libpythia8.a"
    if os.path.exists ( afile ):
        size = os.stat ( afile ).st_size
        if size > 10000000:
            print ( f"[installer.py] {afile} exists. skipping compilation." )
            return True
        else:
            shutil.rmtree ( afile )
    return False

def compilePythia():
    """ finally, compile pythia """
    if checkPythia():
        return
    ver = getVersion()
    ncpus = getNCPUs()
    cmd = f"cd pythia{ver}; ./configure ; make -j {ncpus}"
    print ( f"[installer.py] {cmd}" )
    import subprocess
    ps = subprocess.Popen ( cmd, shell=True, stdout = subprocess.PIPE, \
                            stderr = subprocess.STDOUT )
    for c in iter(lambda: ps.stdout.read(1), b''): 
        sys.stdout.buffer.write(c)
        sys.stdout.buffer.flush()

def installPythia():
    """ fetch tarball, unzip it, compile pythia """
    if checkPythia():
        rmTarball()
        return
    fetch()
    unzip()
    compilePythia()
    rmTarball()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="pythia8 install script" )
    parser.add_argument( '-i', '--install', help='install pythia8',
                         action='store_true') 
    parser.add_argument( '-v', '--version', help='report pythiaversion',
                         action='store_true') 
    args = parser.parse_args()
    if args.version:
        ver = getVersion()
        print ( ver )
        sys.exit()
    installPythia()
