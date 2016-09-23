#!/usr/bin/env python

"""
.. module:: setup
   :synopsis: Setup script for SModelS.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import os
import sys
from setuptools import setup, Extension
sys.path.insert ( 0, "./" )
from smodels.installation import version, authors

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
    ret = [("", [ "README", "COPYING" ])]
    ret.append ( ("smodels/", [ "smodels/version" ]) )
    # ret.append ( ("share", [ "share/shareme" ]) )
    for directory in ["inputFiles/slha/", "inputFiles/lhe/", "lib/nllfast/nllfast-1.2/", "lib/nllfast/nllfast-2.1/", "lib/nllfast/nllfast-3.1/", "lib/pythia6/", "etc/"]:
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


compile()
setup(
    name = "smodels",
    version = version(),
    author = authors(),
    author_email="smodels-developers@lists.oeaw.ac.at ",
    scripts=[ "bin/smodels-config", "runSModelS.py", "smodelsTools.py" ],
    install_requires=[ 'docutils>=0.3', 'numpy', 'scipy>=0.9.0', \
                         'unum', 'argparse', 'pyslha>=3.1.0' ],
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
        "Development Status :: 5 - Production/Stable",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ]
)
