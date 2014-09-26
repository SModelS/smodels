#!/usr/bin/env python

"""
.. module:: testTimeout
   :synopsis: checks if the timeout of pythia6 works
    
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
    
"""

import unittest
import time
import logging.config
from smodels.tools import externalPythia6

class TimeoutTest(unittest.TestCase):
    # use different logging config for the unit tests.
    logging.config.fileConfig( "./logging.conf" )
    logger = logging.getLogger(__name__)
    pythia = externalPythia6.ExternalPythia6()
    pythia.replaceInCfgFile ( {"NEVENTS": 1000, "SQRTS":8000} )

    def testNegative(self):
        t0=time.time()
        self.logger.info ( "test that we dont get timed out" )
        ret=self.pythia.run ( "../inputFiles/slha/T1.slha" )
        lines=ret.split("\n")
        dt=time.time()-t0
        print "dt=",dt
        subprocesses=lines[-12]
        self.assertTrue( dt > .03 )
        self.assertTrue( subprocesses.find("All included")>-1) 
        self.assertTrue( subprocesses.find("1000")>-1) 

    def mestPositive(self):
        t0=time.time()
        self.logger.info ( "test if we do get timed out" )
        print "dt=",time.time()-t0
        self.assertAlmostEqual( 1.0, 1.0 )


if __name__ == "__main__":
    unittest.main()
