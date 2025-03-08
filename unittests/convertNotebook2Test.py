#!/usr/bin/env python3

"""
.. module:: convertNotebook2Test
   :synopsis: Reads the cells from a notebook and store the output from each cell execution.


.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""
import sys
import os
import json
import pickle

def getNotebookData(notebookFile,savePickle=True):

    if not os.path.isfile(notebookFile):
        print('File %s not found' %notebookFile)
        return None

    with open(notebookFile,'r') as f:
        data = json.load(f)

    cells = data['cells']
    outputDict = {}
    for cell in cells:
        cell_type = cell['cell_type']
        if cell_type != 'code':
            continue
        cell_count = cell['execution_count']
        if cell_count is None:
            continue
            
        outputList = [[outStr.replace('\n','') for outStr in c['text']]
                      for c in cell['outputs'] if c['output_type'] == 'stream']
        outputDict[cell_count] = outputList

    if savePickle:
        outputFile = notebookFile
        outputFile = outputFile.replace('.ipynb','_output.pcl')
        pickle.dump(outputDict,open(outputFile,'wb'))
    return outputDict


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser('Convert the notebook output to a pickle file to be used as a unit test')
    ap.add_argument('-f','--filename', help='path to the notebook file')
    args = ap.parse_args()


    getNotebookData(args.filename)
