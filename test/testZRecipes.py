#!/usr/bin/env python3

"""
.. module:: testZRecipes
   :synopsis: Tests that the recipes runs (the Z in the name is,
              so it gets tested last)

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import unittest
import time
import os,glob
import sys
import subprocess
sys.path.append('../')
from smodels.tools.smodelsLogging import logger
from smodels.tools.smodelsLogging import setLogLevel
setLogLevel( "info" )

try:
    from parameterized import parameterized
except:
    logger.error("Install parameterized in order to run notebook tests")
try:
    import pytest
except:
    logger.error("pip install pytest and nbmake for testing notebooks")


nbdir = "../docs/manual/source/recipes"

def listOfNotebooks(notebookDir=nbdir):
    import glob
    notebooks = glob.glob(f"{notebookDir}/*.ipynb")
    notebooks = list(notebooks)
    notebooks.sort()
    notebooks = [os.path.basename(n) for n in notebooks]
    return notebooks

def checkArgParser ():
    import argparse
    try:
        parser = argparse.ArgumentParser ( allow_abbrev=True )
    except TypeError as e:
        print ( "argparser does not take allow_abbrev as argument." )
        print ( "this might be a problem. If it is, consider DOWNgrading argparser!" )

checkArgParser()

class RecipeTest(unittest.TestCase):

    @parameterized.expand(listOfNotebooks())
    def testRunRecipe(self,notebookFile):

        filename = os.path.join(nbdir,notebookFile)
        cmd = f"pytest --nbmake --nbmake-timeout=900 {filename}"
        # print ( "cmd", cmd )
        p = subprocess.Popen([cmd], shell=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        t0 = time.time()
        output, error = p.communicate()
        logger.debug("%s run in %1.2fs" %(notebookFile,time.time()-t0))
        if p.returncode != 0:
            logger.debug("Notebook %s failed with error:\n %s" %(notebookFile,error))
        self.assertEqual(p.returncode,0)



if __name__ == "__main__":
    unittest.main()
