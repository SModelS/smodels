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
    subprocess.call(["jupyter nbconvert --output-dir='./' --execute --to notebook %s" % notebookFile], shell=True)

    outputFile = os.path.basename(notebookFile)
    outputDict = getNotebookData(outputFile, savePickle=False)

    if os.path.isfile(outputFile):
        os.remove(outputFile)

    return outputDict


class Test_Notebook():
    @pytest.mark.parametrize(
        'notebookFile', list(glob.glob('../notebooks-Examples/*.ipynb')))
    def test_notebook(self, notebookFile):

        defaultOutputDict = loadDefaultOutput(notebookFile)
        assert defaultOutputDict is not None
        outputDict = getNewOutput(notebookFile)
        assert outputDict is not None

        default_cells = list(defaultOutputDict.keys())
        new_cells = [icell for icell in outputDict.keys() if icell not in default_cells]
        cells = sorted(default_cells)
        for icell in cells:
            if isinstance(defaultOutputDict[icell], list):
                for iout, out in enumerate(defaultOutputDict[icell]):
                    assert out == outputDict[icell][iout]
            else:
                assert defaultOutputDict[icell] == outputDict[icell]

        assert new_cells == []  # Make sure the notebook has no more cells than the default
