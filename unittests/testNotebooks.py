#!/usr/bin/env python3

"""
.. module:: pytestExampleNotebooks
   :synopsis: Tests the notebooks in the ../notebooks-Examples folder.
              Specific tests can be selected running:
              pytest pytestExampleNotebooks.py -v -k <string>
              where <string> is a string which should be in the notebook filename.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import pickle
import subprocess
from convertNotebook2Test import getNotebookData
import unittest
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
        errorMsg = "Error fetching default output for %s" %notebookFile
        self.assertTrue(defaultOutputDict is not None,msg=errorMsg)

        # Run notebook:
        filename = os.path.join(nbdir,notebookFile)
        p = subprocess.Popen(["jupyter nbconvert --output-dir='./' --execute --to notebook %s" % filename],
                          shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _, error = p.communicate()
        if p.returncode != 0:
            logger.error("Notebook error: %s" %error)
        errorMsg = "Error running notebook %s" %notebookFile
        self.assertEqual(p.returncode,0, msg=errorMsg)

        # Get notebook output
        outputFile = os.path.basename(notebookFile)
        outputDict = getNotebookData(outputFile, savePickle=False)
        errorMsg = "Error executing %s (empty dict)" %notebookFile
        self.assertTrue(outputDict is not None, msg=errorMsg)

        if os.path.isfile(outputFile):
                os.remove(outputFile)

        default_cells = list(defaultOutputDict.keys())
        new_cells = [icell for icell in outputDict.keys() if icell not in default_cells]
        errorMsg = 'New notebook has %i extra cells (cells = %s)' %(len(new_cells),new_cells)
        self.assertTrue(len(new_cells) == 0, msg = errorMsg)  # Make sure the notebook has no more cells than the default
        old_cells = [icell for icell in default_cells if icell not in outputDict.keys()]
        errorMsg = 'Old notebook has %i extra cells (cells = %s)' %(len(old_cells),old_cells)
        self.assertTrue(len(old_cells) == 0, msg = errorMsg)  # Make sure the notebook has no more cells than the default



        cells = sorted(default_cells)
        passed = True
        for icell in cells:
            if isinstance(defaultOutputDict[icell], list):
                for iout, out in enumerate(defaultOutputDict[icell]):
                    errorMsg = 'Error in cell %i and entry %i' %(icell,iout)
                    try:
                        self.assertEqual(out,outputDict[icell][iout], 
                                     msg=errorMsg)
                    except AssertionError as e:
                        passed = False
                        errorMsg += "\n  Error: %s\n" %(e)
                        print(errorMsg)
            else:
                errorMsg = 'Error in cell %i' %icell
                try:
                    self.assertEqual(defaultOutputDict[icell],outputDict[icell], 
                                 msg=errorMsg)
                except AssertionError as e:
                    passed = False
                    errorMsg += "\n  Error: %s\n" %(e)
                    print(errorMsg)

        self.assertTrue(passed,msg="The notebook output does not match the one found in %s" 
                        %(notebookFile.replace('.ipynb', '_output.pcl')))





if __name__ == "__main__":
    unittest.main()
