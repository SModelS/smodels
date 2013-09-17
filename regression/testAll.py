#!/usr/bin/env python

"""
.. module:: testAll
   :synopsis: Runs all test suites.
    
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
    
"""
import unittest


def main():
    """All tests should be added to the list 'testSuites'."""
    
    testSuites = [
        'testXSecComputer',
        'testDescriptions',
        'testMainProgram',
        'testSLHADecomposer',
        'testSLHAUniqueName'
    ]
    
    suite = unittest.TestSuite()
    
    for test in testSuites:
        try:
            # If the module defines a suite() function, call it to get the suite.            
            mod = __import__(test, globals(), locals(), ['suite'])
            suitefn = getattr(mod, 'suite')
            suite.addTest(suitefn())
        except (ImportError, AttributeError):
            # Else, just load all the test cases from the module.
            suite.addTest(unittest.defaultTestLoader.loadTestsFromName(test))
            
    unittest.TextTestRunner().run(suite)


if __name__ == "__main__":
    main()