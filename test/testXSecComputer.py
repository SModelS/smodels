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
        computer = xsecComputer.XSecComputer ( LO, 2000, 6 )
        w = computer.compute(8*TeV, slhafile ).getDictionary()
        # print ( w )
        w8lo= w[(1000021, 1000021)]['8 TeV (LO)'].asNumber( fb )
        self.assertAlmostEqual(w8lo/264., 1., 2 ) 

    def testNLLGlu(self):
        """ test the computation of NLL cross section """
        self.logger.info ( "test NLL xsecs @ 8 TeV" )
        slhafile="./testFiles/slha/simplyGluino.slha"
        computer = xsecComputer.XSecComputer ( NLL, 2000, 6 )
        w = computer.compute( 8*TeV, slhafile ).getDictionary()
        w8nll= w[(1000021, 1000021)]['8 TeV (NLO+NLL)'].asNumber( fb )
        self.assertAlmostEqual(w8nll / 573., 1., 2 )

    def testLOGlu13(self):
        """ test the computation of LO cross section, pythia6 """
        self.logger.info ( "test LO xsecs @ 13 TeV" )
        slhafile="./testFiles/slha/simplyGluino.slha"
        computer = xsecComputer.XSecComputer ( LO, 3000, 6 )
        w = computer.compute( 13*TeV, slhafile ).getDictionary()
        w13lo= w[(1000021, 1000021)]['13 TeV (LO)'].asNumber( fb )
        self.assertAlmostEqual(w13lo / 2237., 1., 1 )

    def testNLLGlu13(self):
        """ test the computation of NLL cross section with pythia6 """
        self.logger.info ( "test NLL xsecs @ 13 TeV" )
        slhafile="./testFiles/slha/simplyGluino.slha"
        computer = xsecComputer.XSecComputer ( NLL, 3000, 6 )
        w = computer.compute( 13*TeV, slhafile ).getDictionary()
        w13nll= w[(1000021, 1000021)]['13 TeV (NLO+NLL)'].asNumber( fb )
        self.assertAlmostEqual(w13nll / 4320. , 1., 2 )
        
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
        xargs.noautocompile = True
        xargs.nevents = 2000
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
        xargs.ssmultipliers = None
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
        self.assertAlmostEqual(lo/264.,1.,2)
        self.assertAlmostEqual(nll/573.6,1.,2)
        #Check 13 TeV xsecs:
        lo = xsecsInfile.getXsecsFor('13 TeV (LO)')[0].value.asNumber(fb)
        nll = xsecsInfile.getXsecsFor('13 TeV (NLL)')[0].value.asNumber(fb)
        self.assertAlmostEqual(lo / 2230., 1., 1 )
        self.assertAlmostEqual(nll / 4308., 1., 1 )

    def testSSMultipliers(self):
        """ test the signal strength multipliers """
        
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
        xargs.noautocompile = True
        xargs.nevents = 5000
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
        xargs.ssmultipliers = { (1000021,1000021): 4. }
        # xargs.ssmultipliers = { 1000021: 2. }
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
        self.assertAlmostEqual(lo/1058.444,1.,1)
        self.assertAlmostEqual(nll/2299.046,1.,1)
        #Check 13 TeV xsecs:
        lo = xsecsInfile.getXsecsFor('13 TeV (LO)')[0].value.asNumber(fb)
        nll = xsecsInfile.getXsecsFor('13 TeV (NLL)')[0].value.asNumber(fb)
        self.assertAlmostEqual(lo/8910.76,1.,1 )
        self.assertAlmostEqual(nll/17215.5906, 1.,1)
        
    def testSSJokers(self):
        """ test the signal strength multipliers, with jokers """
        
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
        xargs.noautocompile = True
        xargs.ncpus = 1
        xargs.nevents = 5000
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
        xargs.ssmultipliers = { ('*100002?','*1000021'): 4. }
        # xargs.ssmultipliers = { 1000021: 2. }
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
        self.assertAlmostEqual(lo/1056.,1.,2)
        self.assertAlmostEqual(nll/2294.,1.,2)
        #Check 13 TeV xsecs:
        lo = xsecsInfile.getXsecsFor('13 TeV (LO)')[0].value.asNumber(fb)
        nll = xsecsInfile.getXsecsFor('13 TeV (NLL)')[0].value.asNumber(fb)
        self.assertAlmostEqual(lo/8910.76,1.,2 )
        self.assertAlmostEqual(nll/17215.5906, 1.,1)

if __name__ == "__main__":
    unittest.main()
