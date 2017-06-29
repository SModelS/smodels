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
import subprocess

def read(fname):
    """
    Simple method to read a file (fname) located in the current folder.

    """
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def listDirectory (dirname):
    """ list all files in directory, skip subdirs """
    if dirname[-1] == "/":
        dirname = dirname[:-1]
    files = os.listdir (dirname)
    ret = []
    for file in files:
        fullname = dirname + "/" + file
        extension = os.path.splitext ( file )[1]
        if os.path.isdir ( fullname ) or \
                extension in [ ".out", ".tgz", ".1" ] or \
                file in [ "Makefile", "README" ]:
            continue
        ret.append ( fullname )
    return ret

def dataFiles ():
    """
    List all config files and binaries

    """
    ret = [("", [ "README.rst", "INSTALLATION.rst", "COPYING" ])]
    ret.append ( ("smodels/", [ "smodels/version" ]) )
    for directory in ["inputFiles/slha/", "inputFiles/lhe/", "smodels/share/",
          "smodels/etc/", "smodels/lib/nllfast/nllfast-1.2/", 
          "smodels/lib/nllfast/nllfast-2.1/", "smodels/lib/nllfast/nllfast-3.1/", 
          "smodels/lib/pythia6/", "smodels/lib/pythia8/" ]:
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
    subprocess.call(["make", "-C", "smodels/lib" ])

# compile() ## not needed anymore as we perform compilation-on-demand now

setup(
    name = "smodels",
    version = version(),
    author = authors(),
    author_email="smodels-developers@lists.oeaw.ac.at ",
    entry_points = {
            'console_scripts': ['smodels-config=smodels.installation:main',
                           'runSModelS.py=smodels.tools.runSModelS:main',
                           'smodelsTools.py=smodels.tools.smodelsTools:main' ]
    },
    install_requires=[ 'docutils>=0.3', 'scipy', 'numpy', 'scipy>=0.9.0', \
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
    include_package_data = True,
    test_suite='test',
    long_description=read('README.rst'),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ]
)
