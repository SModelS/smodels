#!/usr/bin/env python3

"""
.. module:: testCpp
   :synopsis: Tests the C++ interface

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import sys
sys.path.insert(0,"../")
import unittest
import subprocess

class CppTest(unittest.TestCase):
    def compile(self):
        """ compile the C++ interface """
        cmd = "cd ../cpp; make"
        a = subprocess.getoutput ( cmd ).split("\n" )
        ## now check if "done" appears in the last line
        self.assertTrue ( "done" in a[-1] )

    def writeOutput(self):
        """ write the default output """
        f=open("default_cpp.txt","w" )
        cmd = "cd ../cpp; ./run"
        a = subprocess.getoutput ( cmd )
        f.write ( a )
        f.close()

    def runExample(self):
        """ now run the example """
        cmd = "cd ../cpp; ./run"
        a = subprocess.getoutput ( cmd )
        with open("default_cpp.txt","r" ) as f:
            b = f.read()
        self.assertEqual ( a, b )

    def testRun(self):
        self.compile()
        self.runExample()

if __name__ == "__main__":
    unittest.main()
