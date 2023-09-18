#!/usr/bin/env python3

"""
.. module:: testXSecResummino
   :synopsis: Compares the output of XSecResummino with a given value.
    
.. moduleauthor:: Théo Reymermier <theo.reymermier@gmail.com>
    
"""

import sys,os
sys.path.insert(0,"../")
from smodels.tools import xsecComputer, toolBox, xsecResummino
from smodels.tools.xsecComputer import LO, NLL
from smodels.tools.physicsUnits import TeV, fb
from smodels.theory import crossSection
import tempfile
import unittest
import argparse
import logging.config


class XSecTest(unittest.TestCase):
    # use different logging config for the unit tests.
    logging.config.fileConfig( "./logging.conf" )
    from smodels.tools.smodelsLogging import setLogLevel,logger
    
    setLogLevel ( "warn" )

    toolBox.ToolBox().compile() ## make sure the tools are compiled
    
    def testXSecMain(self):
        """ test the main routine for computation of LO and NLL cross section for several sqrts"""
        
        slhafile="./testFiles/resummino/resummino_example.slha"
        f = open(slhafile,'r')
        fdata = f.read()
        fdata = fdata[:fdata.find('XSECTION')]
        f.close()
        fnew = tempfile.mkstemp()
        os.close(fnew[0])
        tmpfile = fnew[1]
        fnew = open(tmpfile,'w')
        fnew.write(fdata)
        fnew.close()        
        self.logger.info ("test NLO xsecs @ 8 and 13 TeV" )
        #Set overall options:
        #Options for cross section calculation:
        xargs = argparse.Namespace()
        xargs.sqrts = [[8.,13.]]
        xargs.ncpus = 1
        xargs.noautocompile = True
        #Compute LO xsecs:
        xargs.query = False
        xargs.conf = 'default'
        xargs.NLL = False
        xargs.NLO = True
        xargs.LOfromSLHA = False
        xargs.keep = False
        xargs.tofile = True
        xargs.alltofile = False
        xargs.pythia6 = True
        xargs.filename = tmpfile
        xargs.colors = False
        xargs.ssmultipliers = None
        xargs.verbosity = "warning"
        #Compute LO cross sections
        xargs.tofile = False
        xargs.alltofile = True
        xsecResummino.main(xargs)
        #Read xsecs:
        xsecsInfile = crossSection.getXsecFromSLHAFile(tmpfile)
        os.remove(tmpfile)
        
        #Check 8 TeV xsecs:
        lo = xsecsInfile.getXsecsFor('8 TeV (LO)')[0].value.asNumber(fb)
        nlo = xsecsInfile.getXsecsFor('8 TeV (NLO)')[0].value.asNumber(fb)
        print('8 TeV is', lo, nlo)
        self.assertAlmostEqual(lo/0.3186,1.,2)
        self.assertAlmostEqual(nlo/0.3691,1.,2)
        #Check 13 TeV xsecs:
        lo = xsecsInfile.getXsecsFor('13 TeV (LO)')[0].value.asNumber(fb)
        nlo = xsecsInfile.getXsecsFor('13 TeV (NLO)')[0].value.asNumber(fb)
        print('13 TeV is', lo, nlo)
        self.assertAlmostEqual(lo / 0.3186, 1., 1 )
        self.assertAlmostEqual(nlo / 0.3691, 1., 1 )


if __name__ == "__main__":
    unittest.main()
