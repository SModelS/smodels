#!/usr/bin/env python

"""
.. module:: testPythia8
   :synopsis: Compares the output of XSecComputer with a given value.
    
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
    
"""

import sys,os
sys.path.insert(0,"../")
from smodels.tools import xsecComputer, toolBox
from smodels.tools.xsecComputer import LO, NLL
from smodels.tools.physicsUnits import TeV, fb
from smodels.theory import crossSection
from unum import Unum
import tempfile
import unittest
import argparse
import logging.config
from math import sqrt

Unum.VALUE_FORMAT = "%.4f"  

Nevents = 50000

def compareXSections(dictA,dictB,nevts,relError = 0.1):

    missingXsecs = set(dictA.keys()).symmetric_difference(set(dictB.keys()))
    commonXsecs = set(dictA.keys()).intersection(set(dictB.keys()))
    totXsecA = sum([x.values()[0].asNumber(fb) for x in dictA.values()])
    totXsecB = sum([x.values()[0].asNumber(fb) for x in dictB.values()])
    totXsec = max(totXsecA,totXsecB)*fb
    mcError = totXsec/sqrt(float(nevts))
    
    diffXsecs = []

    equalXsecs = True
    for xsec in commonXsecs:
        diff = abs(dictA[xsec].values()[0] - dictB[xsec].values()[0])
        if diff > 2.*mcError and (diff/(abs(dictA[xsec].values()[0] + dictB[xsec].values()[0]))).asNumber() > relError:
            print 'Cross-section for %s differ by (%s +- %s) fb' %(str(xsec),str(diff.asNumber(fb)),str(mcError.asNumber(fb)))
            print '   %s   %s' %(dictA[xsec].values()[0],dictB[xsec].values()[0])
            equalXsecs = False

    for xsec in missingXsecs:
        if xsec in dictA:
            if dictA[xsec].values()[0] > 2.*mcError:
                print 'Cross-section for %s missing' %str(xsec)
                print '    %s' %(dictA[xsec].values()[0])
                equalXsecs = False
        else:
            if dictB[xsec].values()[0] > 2.*mcError:
                print 'Cross-section for %s missing' %str(xsec)
                print '    %s' %(dictB[xsec].values()[0])
                equalXsecs =  False
                
    return equalXsecs


class XSecTest(unittest.TestCase):
    # use different logging config for the unit tests.
    logging.config.fileConfig( "./logging.conf" )
    from smodels.tools.smodelsLogging import logger, setLogLevel
    setLogLevel ( "warn" )

    toolBox.ToolBox().compile() ## make sure the tools are compiled

    def testLOGlu(self):
        """ test the computation of LO cross section and compare with pythia6 """
        self.logger.info ("test LO xsecs @ 8 TeV")
        slhafile  = "../inputFiles/slha/gluino_squarks.slha"
        computer6 = xsecComputer.XSecComputer(LO, Nevents, 6)
        computer8 = xsecComputer.XSecComputer(LO, Nevents, 8)
        w6 = computer6.compute(8*TeV, slhafile).getDictionary()
        w8 = computer8.compute(8*TeV, slhafile, pythiacard = './pythia8_to_pythia6.cfg').getDictionary()
        self.assertEqual(compareXSections(w6, w8,Nevents),[]) 
        slhafile  = "../inputFiles/slha/lightEWinos.slha"
        computer6 = xsecComputer.XSecComputer(LO, Nevents, 6)
        computer8 = xsecComputer.XSecComputer(LO, Nevents, 8)
        w6 = computer6.compute(8*TeV, slhafile).getDictionary()
        w8 = computer8.compute(8*TeV, slhafile, pythiacard = './pythia8_to_pythia6.cfg').getDictionary()
        self.assertEqual(compareXSections(w6, w8,Nevents),[]) 
        

if __name__ == "__main__":
    unittest.main()
