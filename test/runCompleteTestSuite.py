#!/usr/bin/env python

"""
.. module:: runCompleteTestSuite
   :synopsis: Runs all test suites.
    
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com> 
    
"""

import unittest
import setPath
import os

def testScript ( filename ):
    """ is filename of the form 'test*py'? """
    return filename[:4]=="test" and filename[-3:]==".py"

def main():
    ## look for test*.py files.
    testSuites = filter(testScript,os.listdir(".") )

    suite = unittest.TestSuite()
    
    for test in testSuites:
        try:
            # If the module defines a suite() function, call it to get the suite.            
            mod = __import__(test[:-3], globals(), locals(), ['suite'])
            suitefn = getattr(mod, 'suite')
            suite.addTest(suitefn())
        except (ImportError, AttributeError):
            # Else, just load all the test cases from the module.
            suite.addTest(unittest.defaultTestLoader.loadTestsFromName(test[:-3]))
            
    unittest.TextTestRunner().run(suite)


if __name__ == "__main__":
    main()
