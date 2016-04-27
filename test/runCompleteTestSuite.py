#!/usr/bin/env python

"""
.. module:: runCompleteTestSuite
   :synopsis: Runs all test suites.
    
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com> 
    
"""

import sys
sys.path.insert(0,"../")
import unittest

def run():
    unittest.TextTestRunner().run( unittest.TestLoader().discover("./") )

def verbose_run():
    alltests = unittest.TestLoader().discover("./") 
    for series in alltests:
        for test in series:
            for t in test:
                print "Running ",t.id()
                t.run()

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser('runs the complete test suite')
    ap.add_argument('-v','--verbose', help='run verbosely',action='store_true')
    args = ap.parse_args()
    if args.verbose:
        verbose_run()
    else:
        run()
