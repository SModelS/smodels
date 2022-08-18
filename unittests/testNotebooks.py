#!/usr/bin/env python3

"""
.. module:: pytestExampleNotebooks
   :synopsis: Tests the notebooks in the ../notebooks-Examples folder.
              Specific tests can be selected running:
              pytest pytestExampleNotebooks.py -v -k <string>
              where <string> is a string which should be in the notebook filename.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import glob
import os
import pickle
import subprocess
from convertNotebook2Test import getNotebookData
import pytest
import unittest
import time
import os,glob
import sys
import subprocess
sys.path.append('../')
from smodels.base.smodelsLogging import logger
from smodels.base.smodelsLogging import setLogLevel
setLogLevel( "info" )

try:
    from parameterized import parameterized
except:
    logger.error("Install parameterized in order to run notebook tests")
try:
    import pytest
except:
    logger.error("pip install pytest and nbmake for testing notebooks")


nbdir = "./notebookTests"


def listOfNotebooks(notebookDir=nbdir):
    import glob
    notebooks = glob.glob(f"{notebookDir}/*.ipynb")
    notebooks = list(notebooks)
    notebooks = [os.path.basename(n) for n in notebooks]
    return notebooks

def loadDefaultOutput(notebookFile):
    # Get default output
    outputFile = notebookFile.replace('.ipynb', '_output.pcl')
    outputFile = os.path.join(nbdir,outputFile)
    if not os.path.isfile(outputFile):
        print('File %s not found. (Try running convertNotebook2Test)' % outputFile)
        return None

    with open(outputFile, 'rb') as f:
        defaultOutputDict = pickle.load(f)
    return defaultOutputDict


class TestNotebook(unittest.TestCase):

    @parameterized.expand(listOfNotebooks())
    def test_notebook(self, notebookFile):

        defaultOutputDict = loadDefaultOutput(notebookFile)
        self.assertTrue(defaultOutputDict is not None)

        # Run notebook:
        filename = os.path.join(nbdir,notebookFile)
        p = subprocess.Popen(["jupyter nbconvert --output-dir='./' --execute --to notebook %s" % filename],
                          shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = p.communicate()
        if p.returncode != 0:
            logger.error("Notebook error: %s" %error)
        self.assertEqual(p.returncode,0)

        # Get notebook output
        outputFile = os.path.basename(notebookFile)
        outputDict = getNotebookData(outputFile, savePickle=False)
        self.assertTrue(outputDict is not None)

        if os.path.isfile(outputFile):
                os.remove(outputFile)

        default_cells = list(defaultOutputDict.keys())
        new_cells = [icell for icell in outputDict.keys() if icell not in default_cells]
        self.assertTrue(len(new_cells) == 0)  # Make sure the notebook has no more cells than the default

        cells = sorted(default_cells)
        for icell in cells:
            # print(icell)
            # print(defaultOutputDict[icell])
            # print(outputDict[icell],'\n')
            if isinstance(defaultOutputDict[icell], list):
                for iout, out in enumerate(defaultOutputDict[icell]):
                    self.assertEqual(out,outputDict[icell][iout])
            else:
                self.assertEqual(defaultOutputDict[icell],outputDict[icell])






if __name__ == "__main__":
    unittest.main()
