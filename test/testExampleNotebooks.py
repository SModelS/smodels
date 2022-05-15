#!/usr/bin/env python3

"""
.. module:: convertNotationTest
   :synopsis: Tests the theory.element.Element class

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""
import sys,glob
sys.path.insert(0,"../")
import os
import unittest
import pickle
import subprocess
from convertNotebook2Test import getNotebookData

def loadDefaultOutput(notebookFile):
    # Get default output
    outputFile = os.path.basename(notebookFile)
    outputFile = outputFile.replace('.ipynb','_output.pcl')
    if not os.path.isfile(outputFile):
        print('File %s not found. (Try running convertNotebook2Test)' %outputFile)
        return None

    with open(outputFile,'rb') as f:
        defaultOutputDict = pickle.load(f)
    return defaultOutputDict

def getNewOutput(notebookFile):
    # Run notebook
    r = subprocess.call(["jupyter nbconvert --output-dir='./' --execute --to notebook %s" %notebookFile],shell=True)

    outputFile = os.path.basename(notebookFile)
    outputDict = getNotebookData(outputFile,savePickle=False)

    if os.path.isfile(outputFile):
        os.remove(outputFile)

    return outputDict


class notebookTests(unittest.TestCase):

    def testNotebooks(self):

        for notebookFile in glob.glob('../notebooks-Examples/*.ipynb'):
            defaultOutputDict = loadDefaultOutput(notebookFile)
            self.assertFalse(defaultOutputDict is None)
            outputDict = getNewOutput(notebookFile)
            self.assertFalse(outputDict is None)

            print('\nChecking %s' %os.path.basename(notebookFile))
            for icell in sorted(defaultOutputDict.keys()):
                print('\t checking cell %i' %icell)
                self.assertEqual(defaultOutputDict[icell],outputDict[icell])

if __name__ == "__main__":
    unittest.main()
