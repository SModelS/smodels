#!/usr/bin/env python3

"""
.. module:: testCpp
   :synopsis: Tests the C++ interface

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import sys
sys.path.insert(0,"../")
from smodels.tools.smodelsLogging import setLogLevel
setLogLevel ( "error" )

import unittest
if sys.version[0]=="2":
    import commands as CMD
else:
    import subprocess as CMD

class CppTest(unittest.TestCase):
    def compile(self):
        """ compile the C++ interface """
        cmd = "cd ../cpp; make"
        a = CMD.getoutput ( cmd ).split("\n" )
        ## now check if "done" appears in the last line
        self.assertTrue ( "done" in a[-1] )

    def writeOutput(self):
        """ write the default output. will *define* the unit test. """
        f=open("default_cpp.txt","w" )
        cmd = "cd ../cpp; ./run"
        a = CMD.getoutput(cmd)
        f.write ( a )
        f.close()

    def writeIni(self, dbpath ):
        """ write the ini file """
        templf = open ( "../cpp/template.ini", "rt" )
        templines = templf.readlines()
        templf.close()
        if "database" in dbpath and not "https" in dbpath:
                dbpath = "../test/"+dbpath
        parmf = open ( "../cpp/parameters.ini", "wt" )
        for line in templines:
            parmf.write ( line.replace ( "@@DBPATH@@", dbpath ) )
        parmf.close()

    def runExample(self):
        """ now run the example """
        from databaseLoader import dbpath
        self.writeIni ( dbpath )
        cmd = "cd ../cpp; ./run"
        l = CMD.getoutput(cmd)
        l = l[l.find('Input status'):]
        la = l.split("\n")
        a = [ x.strip() for x in la  if x.strip() and x.strip()[0] != '#' 
             and not 'WARNING' in x and not 'INFO' in x]        
        with open("default_cpp.txt","r" ) as f:
            l = f.read()
            l = l[l.find('Input status'):]            
            lb = l.split("\n")
            b = [ x.strip() for x in lb  if x.strip() and x.strip()[0] != '#' 
                 and not 'WARNING' in x and not 'INFO' in x]
        if len(a) != len(b):
            print ( "test failed. writing output to debug.txt" )
            f=open("debug.txt","w")
            for i in b:
                f.write ( i+ "\n" )
            f.close()
        self.assertEqual(len(a), len(b))
        for x,y in zip(a,b):
            if x == y:
                continue
            xvals = [v.strip() for v in x.split()]
            yvals = [v.strip() for v in y.split()]
            self.assertEqual(len(xvals), len(yvals), "Lines:\n %s and \n %s \n differ" %(x,y))
            for i,xv in enumerate(xvals):
                yv = yvals[i]
                if "=======" in xv and "=======" in yv:
                    continue
                try:
                    yv = eval(yv)
                    xv = eval(xv)
                except:
                    pass
                if isinstance(yv,float) and isinstance(xv,float):
                    self.assertAlmostEqual(xv, yv, 3)
                else:
                    self.assertEqual(xv, yv)
        ## reset parameter.ini file
        self.writeIni ( "unittest" )

    def testRun(self):
        self.compile()
        #self.writeOutput()
        self.runExample()

if __name__ == "__main__":
    unittest.main()
