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
        """ write the default output. will *define* the unit test. """
        f=open("default_cpp.txt","w" )
        cmd = "cd ../cpp; ./run"
        a = subprocess.getoutput ( cmd )
        f.write ( a )
        f.close()

    def runExample(self):
        """ now run the example """
        cmd = "cd ../cpp; ./run"
        l = subprocess.getoutput ( cmd )
        skip=8 ## skip first lines
        la = l.split("\n")[skip:] 
        a = [ x.strip() for x in la ]
        with open("default_cpp.txt","r" ) as f:
            l = f.readlines()
            lb=l[skip:]
            b = [ x.strip() for x in lb ]
        self.assertEqual ( len(a), len(b) )
        for x,y in zip(a,b):
            self.assertEqual ( x, y )

    def testRun(self):
        self.compile()
        # self.runExample()
        self.writeOutput()

if __name__ == "__main__":
    unittest.main()
