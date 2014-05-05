#!/usr/bin/env python

"""
.. module:: callGraph
   :synopsis: Call pycallgraph to create call graphs of simpleExample for
   different input files.
   
   Usage: has to be invoked in the SModelS root directory:
   
   >>> python -m tools.callGraph
   >>> python -m tools.callGraph --max-depth 3

.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>

"""

from __future__ import print_function
import argparse
import setPath
from smodels import simpleExample
from os.path import os

inputFiles = ['inputFiles/slha/andrePT4.slha',
              'inputFiles/slha/T1.slha']


def main():
    """
    Executes pycallgraph to create call graphs of simpleExample.
    
    """
    synopsis = "Call pycallgraph to create call graphs of simpleExample."
    desc = ("Maximum stack depth to trace. Any calls made past this stack "
            "depth are not included in the trace.")

    parser = argparse.ArgumentParser(description=synopsis)
    parser.add_argument('--max-depth', type=int, choices=[2, 3, 4, 5],
                        help=desc)
    args = parser.parse_args()

    depth = ""
    if args.max_depth:
        depth = '--max-depth ' + str(args.max_depth) + ' '

    arguments = ('--include "theory.*" --include "experiment.*" '
                 '--include "tools.*" ' + depth + 'graphviz -- '
                 './simpleExample.py')

    for inputFile in inputFiles:
        simpleExample.slhafile = inputFile

        print("calling pycallgraph:")
        print("    input file:", inputFile)
        print("     arguments:", str(arguments))

        os.system('pycallgraph ' + arguments)

        filename = os.path.splitext(os.path.basename(inputFile))[0] + ".png"

        print("renaming output file to", filename)

        cwd = os.getcwd()
        oldFile = os.path.join(cwd, 'pycallgraph.png')
        newFile = os.path.join(cwd, filename)
        print(oldFile)
        print(newFile)
        if os.path.isfile(oldFile):
            os.rename(oldFile, newFile)
        else:
            print("ERROR:", oldFile, "does not exist!")


if __name__ == '__main__':
    main()
