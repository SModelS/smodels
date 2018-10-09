#!/usr/bin/env python3

"""
.. module:: testCpp
   :synopsis: Tests the C++ interface

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import sys
sys.path.insert(0,"../")
from smodels.tools.smodelsLogging import setLogLevel

import unittest
if sys.version[0]=="2":
    import commands as CMD
else:
    import subprocess as CMD

class CppTest(unittest.TestCase):
    skip=8 ## skip first lines

    def compile(self):
        """ compile the C++ interface """
        cmd = "cd ../cpp; make"
        a = CMD.getoutput ( cmd ).split("\n" )
        ## now check if "done" appears in the last line
        self.assertTrue ( "done" in a[-1] )

    def writeOutput(self):
        """ write the default output. will *define* the unit test. """
        setLogLevel ( "error" )
        f=open("default_cpp.txt","w" )
        cmd = "cd ../cpp; ./run"
        a = CMD.getoutput ( cmd )
        b = a.split("\n")
        for ctr,i in enumerate(b):
            if "INFO" in i or "WARNING" in i or "smodels.cpp" in i:
                continue
            ctr+=1
            if ctr < self.skip:
                continue
            f.write ( i+ "\n" )
        f.close()

    def runExample(self):
        """ now run the example """
        setLogLevel ( "error" )
        cmd = "cd ../cpp; ./run"
        l = CMD.getoutput ( cmd )
        la = l.split("\n")[self.skip:] 
        a = [ x.strip() for x in la ]
        with open("default_cpp.txt","r" ) as f:
            l = f.readlines()
            lb=l ## [self.skip:]
            b = [ x.strip() for x in lb ]
        if len(a) != len(b):
            print ( "test failed. writing output to debug.txt" )
            f=open("debug.txt","w")
            for i in b:
                if "INFO" in i or "WARNING" in i or "smodels.cpp" in i:
                    continue
                f.write ( i+ "\n" )
            f.close()
        self.assertEqual ( len(a), len(b) )
        for x,y in zip(a,b):
            self.assertEqual ( x, y )

    def testRun(self):
        self.compile()
        # self.writeOutput()
        self.runExample()

if __name__ == "__main__":
    unittest.main()
