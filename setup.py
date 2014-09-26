#!/usr/bin/env python

import os
from setuptools import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    """
    Simple method to read a file (fname) located in the current folder.
    """
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

# nllfast={"1.2": "nllfast_7TeV", "2.1": "nllfast_8TeV", "3.0":"nllfast_13TeV",
#          "4.01dcpl": "nllfast_14TeV", "5.01dcpl": "nllfast_33TeV" }

def dataFiles ():
    """ obtain all data files dynamically """
    ret = [("etc/", ["etc/logging.conf", "etc/pythia.card"])]
    from smodels.tools import toolBox
    box = toolBox.ToolBox()
    print box.listOfTools()
    # for (k,v) in nllfast.items():
    #    pth="tools/external/nllfast/nllfast-%s/" % k
    #    ret.append( ( pth, [  pth + v ] ) )
    #    for fle in os.listdir ( pth ):
    #        ## if fle in [ "README", "Makefile" ]: continue
    #        if fle[-5:]!=".grid": continue
    #        ## print fle
    #        ret.append( ( pth, [  pth + fle ] ) )
    return ret

setup(
    name = "smodels",
    version = "1.0",
    author = ("Sabine Kraml, Suchita Kulkarni, Ursula Laa, Andre Lessa, "
              "Wolfgang Magerl, Doris Proschofsky, Wolfgang Waltenberger"),
    author_email = "smodels-developers@lists.oeaw.ac.at ",
    scripts = ["bin/smodels-config",
               "tools/external/nllfast/nllfast-4.01dcpl/nllfast_14TeV"],
    install_requires = ['docutils>=0.3', 'numpy', 'scipy>=0.9.0', 'unum'],
    data_files = dataFiles() ,
    description = ("A tool for interpreting simplified-model results from the "
                   "LHC"),
    license = "GPLv3",
    # use_2to3 = True,
    keywords = ("simplified models LHC BSM theories interpretation "
                "supersymmetry universal extra dimensions"),
    url = "http://smodels.hephy.at/",
    packages = ['smodels',
                'smodels.theory',
                'smodels.tools',
                'smodels.experiment'],
    test_suite = 'test',
    long_description = read('README'),
    classifiers = [
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ],
)
