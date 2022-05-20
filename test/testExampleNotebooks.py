#!/usr/bin/env python3

"""
.. module:: convertNotationTest
   :synopsis: Tests the theory.element.Element class

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""
import glob
import os
import unittest
import time
import pickle
import subprocess
from convertNotebook2Test import getNotebookData
import progressbar


def loadDefaultOutput(notebookFile):
    # Get default output
    outputFile = os.path.basename(notebookFile)
    outputFile = outputFile.replace('.ipynb', '_output.pcl')
    if not os.path.isfile(outputFile):
        print('File %s not found. (Try running convertNotebook2Test)' % outputFile)
        return None

    with open(outputFile, 'rb') as f:
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

            # if 'decomp' in notebookFile:
                # continue
            # if not 'split' in notebookFile:
                # continue

            defaultOutputDict = loadDefaultOutput(notebookFile)
            self.assertFalse(defaultOutputDict is None)
            outputDict = getNewOutput(notebookFile)
            self.assertFalse(outputDict is None)

            print('\nChecking %s' % os.path.basename(notebookFile))
            widgets = [progressbar.Percentage(), progressbar.Bar(),
                       progressbar.Counter()]
            cells = list(defaultOutputDict.keys()) + list(outputDict.keys())
            cells = sorted(list(set(cells)))
            bar = progressbar.ProgressBar(widgets=widgets,
                                          maxval=len(cells)).start()

            for icell in cells:
                bar.update(icell)
                time.sleep(0.2)
                # print('\t checking cell %i' %icell)
                # print(defaultOutputDict[icell])
                try:
                    self.assertEqual(defaultOutputDict[icell],outputDict[icell])
                except AssertionError:
                    print('----- CELL %i ------' % icell)
                    print('default:\n %s \n new:\n %s\n' %(defaultOutputDict[icell],outputDict[icell]))
                    raise AssertionError()
            bar.finish()

if __name__ == "__main__":
    unittest.main()
