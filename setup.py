#!/usr/bin/env python

import os
from setuptools import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "smodels",
    version = "1.0.0",
    author = "Sabine Kraml, Suchita Kulkarni, Ursula Laa, Andre Lessa, Wolfgang Magerl, Doris Proschofsky, Wolfgang Waltenberger",
    author_email = "smodels-developers@lists.oeaw.ac.at ",
    scripts = [ "bin/smodels-config" ],
    install_requires = ['docutils>=0.3', 'numpy', 'scipy', 'unum' ],
    data_files = [ ( "etc/", [ "smodels/etc/logging.conf" ] ) ],
    description = ("A tool for interpreting simplified-model results from the LHC"),
    license = "GPLv3",
    # use_2to3 = True,
    keywords = "simplified models LHC BSM theories interpretation supersymmetry universal extra dimensions",
    url = "http://smodels.hephy.at/",
    # packages=['experiment', 'theory', 'tools' ],
    packages=[ 'smodels', 'smodels.theory', 'smodels.tools', 'smodels.experiment' ],
    test_suite = 'smodels.tests',
    long_description=read('README'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ],
)
