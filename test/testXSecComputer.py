#!/usr/bin/env python3

"""
.. module:: testXSecComputer
   :synopsis: Compares the output of XSecComputer with a given value.
    
.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
    
"""

import sys,os
sys.path.insert(0,"../")
from smodels.tools import xsecComputer, toolBox
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

    def testLOGlu(self):
        """ test the computation of LO cross section, pythia6 """
        self.logger.info ( "test LO xsecs @ 8 TeV" )
        slhafile= "./testFiles/slha/simplyGluino.slha"
        computer = xsecComputer.XSecComputer ( LO, 10000, 6 )
        w = computer.compute(8*TeV, slhafile ).getDictionary()
        # print ( w )
        w8lo= w[(1000021, 1000021)]['8 TeV (LO)'].asNumber( fb )
        self.assertAlmostEqual(w8lo, 262.453, 2 ) 

    def testNLLGlu(self):
        """ test the computation of NLL cross section """
        self.logger.info ( "test NLL xsecs @ 8 TeV" )
        slhafile="./testFiles/slha/simplyGluino.slha"
        computer = xsecComputer.XSecComputer ( NLL, 10000, 6 )
        w = computer.compute( 8*TeV, slhafile ).getDictionary()
        w8nll= w[(1000021, 1000021)]['8 TeV (NLO+NLL)'].asNumber( fb )
        self.assertAlmostEqual(w8nll, 570.0737, 2 )

    def testLOGlu13(self):
        """ test the computation of LO cross section, pythia6 """
        self.logger.info ( "test LO xsecs @ 13 TeV" )
        slhafile="./testFiles/slha/simplyGluino.slha"
        computer = xsecComputer.XSecComputer ( LO, 10000, 6 )
        w = computer.compute( 13*TeV, slhafile ).getDictionary()
        w13lo= w[(1000021, 1000021)]['13 TeV (LO)'].asNumber( fb )
        self.assertAlmostEqual(w13lo, 2234.6808, 2 )

    def testNLLGlu13(self):
        """ test the computation of NLL cross section with pythia6 """
        self.logger.info ( "test NLL xsecs @ 13 TeV" )
        slhafile="./testFiles/slha/simplyGluino.slha"
        computer = xsecComputer.XSecComputer ( NLL, 10000, 6 )
        w = computer.compute( 13*TeV, slhafile ).getDictionary()
        w13nll= w[(1000021, 1000021)]['13 TeV (NLO+NLL)'].asNumber( fb )
        self.assertAlmostEqual(w13nll, 4317.4, 2 )
        
    def testXSecMain(self):
        """ test the main routine for computation of LO and NLL cross section for several sqrts"""
        
        slhafile="./testFiles/slha/simplyGluino.slha"
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
        self.logger.info ("test NLL xsecs @ 8 and 13 TeV" )
        #Set overall options:
        #Options for cross section calculation:
        xargs = argparse.Namespace()
        xargs.sqrts = [[8.,13.]]
        xargs.ncpus = 1
        xargs.nevents = 10000
        #Compute LO xsecs:
        xargs.query = False
        xargs.NLL = False
        xargs.NLO = False
        xargs.LOfromSLHA = False
        xargs.keep = False
        xargs.tofile = True
        xargs.alltofile = False
        xargs.pythia6 = True
        xargs.filename = tmpfile
        xargs.colors = False
        xargs.verbosity = "warning"
        #Compute LO cross sections
        xsecComputer.main(xargs)
        #Compute NLL cross sections
        xargs.NLL = True
        xargs.LOfromSLHA = True
        xsecComputer.main(xargs)
        #Read xsecs:
        xsecsInfile = crossSection.getXsecFromSLHAFile(tmpfile)
        os.remove(tmpfile)
        
        #Check 8 TeV xsecs:
        lo = xsecsInfile.getXsecsFor('8 TeV (LO)')[0].value.asNumber(fb)
        nll = xsecsInfile.getXsecsFor('8 TeV (NLL)')[0].value.asNumber(fb)
        self.assertAlmostEqual(lo,262.4527,2)
        self.assertAlmostEqual(nll,570.0736,2)
        #Check 13 TeV xsecs:
        lo = xsecsInfile.getXsecsFor('13 TeV (LO)')[0].value.asNumber(fb)
        nll = xsecsInfile.getXsecsFor('13 TeV (NLL)')[0].value.asNumber(fb)
        self.assertAlmostEqual(lo,2234.6808,2 )
        self.assertAlmostEqual(nll, 4317.40338, 2 )

        
        
#         self.assertAlmostEqual(w13nll, 4341.659372400062)
        

if __name__ == "__main__":
    unittest.main()
