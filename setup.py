#!/usr/bin/env python

"""
.. module:: setup
   :synopsis: Setup script for SModelS.
   
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import os
from setuptools import setup, Extension


def read(fname):
    """
    Simple method to read a file (fname) located in the current folder.
    
    """
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def listDirectory (dirname):
    if dirname[-1] == "/":
        dirname = dirname[:-1]
    files = os.listdir (dirname)
    ret = []
    for file in files:
        ret.append (dirname + "/" + file)
    return ret

def dataFiles ():
    """
    List all config files and binaries

    """
    ret = [("", [ "BANNER", "README", "COPYING" ])]
    ret.append ( ( "smodels/", [ "smodels/version" ] ) )

    for directory in ["inputFiles/slha/", "inputFiles/lhe/", "lib/nllfast/nllfast-1.2/", "lib/nllfast/nllfast-2.1/", "lib/pythia6/", "etc"]:
        ret.append ((directory, listDirectory (directory)))

    return ret


def compile():
    """
    Compile external tools by calling make

    """
    import sys
    if len(sys.argv) < 2:
        return
    needs_build = False
    for i in sys.argv[1:]:
        if i in ["build", "build_ext", "build_clib", "install", "install_lib", "bdist", "bdist_rpm", "bdist_dumb", "bdist_wininst", "bdist_wheel", "develop"]:
            needs_build = True
    if not needs_build:
        return
    import subprocess
    subprocess.call(["make", "-C", "lib" ])


def version():
    with open("smodels/version") as f: return f.readline()

compile()
setup(
    name = "smodels",
    version = version(),
    author = ("Sabine Kraml, Suchita Kulkarni, Ursula Laa, Andre Lessa, "
              "Veronika Magerl, Wolfgang Magerl, Doris Proschofsky, "
              "Michael Traub, Wolfgang Waltenberger"),
    author_email="smodels-developers@lists.oeaw.ac.at ",
    scripts=[ "bin/smodels-config", "runSModelS.py" ],
    install_requires=[ 'docutils>=0.3', 'numpy', 'scipy>=0.9.0', \
                         'unum', 'argparse'],
    data_files=dataFiles() ,
    description=("A tool for interpreting simplified-model results from the "
                   "LHC"),
    license="GPLv3",
    # use_2to3 = True,
    keywords=("simplified models LHC BSM theories interpretation "
                "supersymmetry UEDs"),
    url="http://smodels.hephy.at/",
    packages=['smodels',
                'smodels.theory',
                'smodels.tools',
                'smodels.experiment'],
    test_suite='test',
    long_description=read('README'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ]
)
