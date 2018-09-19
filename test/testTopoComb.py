#!/usr/bin/env python3
 
"""
.. module:: testTopoComb
   :synopsis: Tests the combinatin of topologies, for combined results

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
 
"""
 
import sys,os
sys.path.insert(0,"../")
import unittest
from databaseLoader import database
 
# from smodels.tools.smodelsLogging import logger, setLogLevel
from smodels.theory.theoryPrediction import theoryPredictionsFor
from smodels.theory.decomposer import decompose
from smodels.tools.physicsUnits import fb
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.theory.model import Model
 
class TopoCombTest(unittest.TestCase):
    def createSLHAFile(self,case="T1tttt"):
        """ create slha file. case defines which case to create
        (T1,T5,mixed ) """

        ratios = { "T1": 0., "T5": 1. }
        if case == "T1":
            ratios = { "T1": 1., "T5": .0 }
        if case == "mixed":
            ratios = { "T1": .5, "T5": .5 }

        fname = "%s.slha" % case
        f=open( fname,"w")
        f.write ( """BLOCK MASS  # Mass Spectrum
   1000006     120.           # ~t_1
   1000021     1200.          # ~g
   1000022     100.           # ~chi_10
DECAY   1000022     0.00000000E+00   # neutralino1 decays
DECAY   1000006     1.00000000E+00   # stop1 decays
     1.00000000E+00    2     1000022         4   # BR(~st1 -> ~chi_10 c)
DECAY   1000021     1.00000000E+00   # gluino decays
      %.2f              3     1000022        -6         6   # BR(~gl -> N1 tbar t)               
      %.2f              2     1000006        -6   # BR(~g -> ~t_1  tb)
XSECTION  1.30E+04  2212 2212 2 1000021 1000021 # 10000 events, [pb], pythia8 for LO
  0  2  0  0  0  0    0.07           SModelSv1.1.3
""" % ( ratios["T1"], ratios["T5"] ) )
        f.close()
        return fname


    def testCombinedResult(self):
        predXSecs,rvalues={},{}
        for case in [ "T1", "T5", "mixed" ]:
            filename = self.createSLHAFile( case=case )
            model = Model(BSMList,SMList,filename)
            model.updateParticles()            
            deco = decompose ( model )    
                    
            expRes = database.getExpResults( analysisIDs = [ "CMS-SUS-16-050-agg" ] )[0]
            # print ( "Experimental result: %s" % expRes )
            tp = theoryPredictionsFor ( expRes, deco, useBestDataset=False, combinedResults=True )
            for t in tp:
                predXSecs[case]=t.xsection.value
                rvalues[case]=t.getRValue(expected=True)
            if True:
                os.unlink ( filename )
        ## first test: the theory prediction of the mixed scenario should be 25% of the sum
        ## 25%, because the total cross section is a fixed quantity, and half of the mixed scenario
        ## goes into asymmetric branches which we miss out on.
        self.assertAlmostEqual ( (predXSecs["T1"]+predXSecs["T5"]).asNumber(fb), (4*predXSecs["mixed"]).asNumber(fb), 2 )

        ## second test: the r value of the mixed scenario * 2 must be between the r values
        ## of the pure scenarios. The factor of two comes from the fact, that we loose 50% 
        ## to asymmetric branches
        self.assertTrue ( rvalues["T5"] < 2*rvalues["mixed"] < rvalues["T1"] )

 
if __name__ == "__main__":
    unittest.main()
