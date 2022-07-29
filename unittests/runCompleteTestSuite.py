#!/usr/bin/env python3

"""
.. module:: runCompleteTestSuite
   :synopsis: Runs all test suites.

.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from __future__ import print_function
import sys, subprocess
sys.path.insert(0,"../")
from smodels.tools.colors import colors
colors.on = True

v=sys.version_info
if v[0] > 2 or ( v[0]==2 and v[1] > 6 ):
    import unittest
if v[0]==2 and v[1] < 7 and v[1] > 3:
    try:
        import unittest2 as unittest
    except ImportError as e:
        print ( "Error: python v",sys.version,"needs unittest2. Please install." )
        sys.exit()
from smodels.tools.smodelsLogging import setLogLevel
setLogLevel ( "error" )

def run():
    unittest.TextTestRunner().run( unittest.TestLoader().discover("./") )

def verbose_run( flter ):
    alltests = unittest.TestLoader().discover("./")
    n_tests, n_failed = 0, 0
    for series in alltests:
        for test in series:
            if type(test)!=unittest.suite.TestSuite:
                print ( "%sError: could not import ``%s'' {%s}" % \
                        ( colors.error, test, type(test) ) )
                print ( test._exception, colors.reset )
                continue
            for t in test:
                if flter and (not flter in str(t)):
                    continue
                n_tests += 1
                print ( "[#%3d] %s ... " % ( n_tests, t.id() ), end="" )
                sys.stdout.flush()
                try:
                    a=t.debug()
                except Exception as e:
                    n_failed += 1
                    print ( "%s FAILED: %s,%s%s" % \
                            ( colors.error, type(e), str(e), colors.reset ) )
                    continue
                print ( "%sok%s" % ( colors.info, colors.reset ) )

                #a=t.run() ## python3
                # print ( "a=",a )
    print( "[runCompleteTestSuite] %d/%d tests failed." % ( n_failed, n_tests ))

def parallel_run ( verbose ):
    if verbose:
        print ("[runCompleteTestSuite] verbose run not implemented "
               "for parallel version" )
        return
    try:
        from concurrencytest import ConcurrentTestSuite, fork_for_tests
    except ImportError as e:
        print ( "Need to install the module concurrencytest." )
        print ( "pip install --user concurrencytest" )
        return
    from smodels.tools import runtime
    suite = unittest.TestLoader().discover("./")
    ncpus = runtime.nCPUs()
    ## "shuffle" the tests, so that the heavy tests get distributed
    ## more evenly among threads (didnt help, so I commented it out)
    #suite._tests = [ item for sublist in [ suite._tests[x::ncpus] \
    #    for x in range(ncpus) ] for item in sublist ]
    concurrent_suite = ConcurrentTestSuite(suite, fork_for_tests( ncpus ))
    runner = unittest.TextTestRunner()
    runner.run(concurrent_suite)

def cleanDatabase ():
    """ remove database pickle files """
    cmd = "rm -r database/*.pcl database/*TeV/*/*/.*pcl"
    o =subprocess.getoutput ( cmd )

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser('runs the complete test suite')
    ap.add_argument('-c','--clean_database', help='remove database pickle files',
                    action='store_true')
    ap.add_argument('-v','--verbose', help='run verbosely',action='store_true')
    ap.add_argument('-f','--filter', help='run only tests that have <FILTER> in name. Works only with verbose and not parallel. case sensitive.',type=str,default=None)
    ap.add_argument('-p','--parallel', help='run in parallel',action='store_true')
    args = ap.parse_args()
    if args.clean_database:
        cleanDatabase ()
    if args.parallel:
        parallel_run ( args.verbose )
        sys.exit()
    if args.verbose:
        verbose_run( args.filter )
        sys.exit()
    run()
