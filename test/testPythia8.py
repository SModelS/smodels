#!/usr/bin/env python3

"""
.. module:: testPythia8
   :synopsis: Compares the output of XSecComputer with a given value.
    
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
    
"""

from __future__ import print_function
import sys
sys.path.insert(0,"../")
from smodels.tools import xsecComputer, toolBox
from smodels.tools.xsecComputer import LO
from smodels.tools.physicsUnits import TeV, fb
from unum import Unum
import unittest
from math import sqrt

Unum.VALUE_FORMAT = "%.4f"  

Nevents = 10000

def compareXSections(dictA,dictB,nevts,relError = 0.1):

    missingXsecs = set(dictA.keys()).symmetric_difference(set(dictB.keys()))
    commonXsecs = set(dictA.keys()).intersection(set(dictB.keys()))
    totXsecA = sum([list(x.values())[0].asNumber(fb) for x in dictA.values()])
    totXsecB = sum([list(x.values())[0].asNumber(fb) for x in dictB.values()])
    totXsec = max(totXsecA,totXsecB)*fb
    mcError = totXsec/sqrt(float(nevts))
    
    equalXsecs = True
    for xsec in commonXsecs:
        diff = abs(list(dictA[xsec].values())[0] - list(dictB[xsec].values())[0])
        if diff > 2.*mcError and (diff/(abs(list(dictA[xsec].values())[0] + list(dictB[xsec].values())[0]))).asNumber() > relError:
            print ( 'Cross-section for %s differ by (%s +- %s) fb' %(str(xsec),str(diff.asNumber(fb)),str(mcError.asNumber(fb))) )
            print ( '   %s   %s' %(list(dictA[xsec].values())[0],list(dictB[xsec].values())[0]) )
            equalXsecs = False

    for xsec in missingXsecs:
        if xsec in dictA:
            if list(dictA[xsec].values())[0] > 2.*mcError:
                print ( 'Cross-section for %s missing' %str(xsec) )
                print ( '    %s' %(list(dictA[xsec].values())[0]) )
                equalXsecs = False
        else:
            if list(dictB[xsec].values())[0] > 2.*mcError:
                print ( 'Cross-section for %s missing' %str(xsec) )
                print ( '    %s' %(list(dictB[xsec].values())[0]) )
                equalXsecs =  False
                
    return equalXsecs


class XSecTest(unittest.TestCase):
    # use different logging config for the unit tests.

    toolBox.ToolBox().compile() ## make sure the tools are compiled

    def testLOGlu(self):
        """ test the computation of LO cross section and compare with pythia6 """

        slhafile  = "./testFiles/slha/gluino_squarks.slha"
        computer6 = xsecComputer.XSecComputer(LO, Nevents, 6)
        computer8 = xsecComputer.XSecComputer(LO, Nevents, 8)
        w6 = computer6.compute(8*TeV, slhafile).getDictionary()
        w8 = computer8.compute(8*TeV, slhafile, pythiacard = './pythia8_to_pythia6.cfg').getDictionary()
        self.assertEqual(compareXSections(w6, w8,Nevents),True) 
        slhafile  = "./testFiles/slha/lightEWinos.slha"
        computer6 = xsecComputer.XSecComputer(LO, Nevents, 6)
        computer8 = xsecComputer.XSecComputer(LO, Nevents, 8)
        w6 = computer6.compute(8*TeV, slhafile).getDictionary()
        w8 = computer8.compute(8*TeV, slhafile, pythiacard = './pythia8_to_pythia6.cfg').getDictionary()
        self.assertEqual(compareXSections(w6, w8,Nevents),True) 
        

if __name__ == "__main__":
    unittest.main()
