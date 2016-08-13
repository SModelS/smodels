#!/usr/bin/env python

"""
.. module:: runCompleteTestSuite
   :synopsis: Runs all test suites.
    
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com> 
    
"""

import sys
sys.path.insert(0,"../")
v=sys.version_info
if v[0] > 2 or ( v[0]==2 and v[1] > 6 ):
    import unittest
if v[0]==2 and v[1] < 7 and v[1] > 3:
    try:
        import unittest2 as unittest
    except ImportError,e:
        print "Error: python v",sys.version,"needs unittest2. Please install."
        sys.exit()
from logging import ERROR
from smodels.tools.printer import logger as pl
pl.setLevel ( ERROR )


def run():
    unittest.TextTestRunner().run( unittest.TestLoader().discover("./") )

def verbose_run():
    alltests = unittest.TestLoader().discover("./") 
    for series in alltests:
        for test in series:
            for t in test:
                print "[runCompleteTestSuite] now run",t.id()
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
