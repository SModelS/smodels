#!/usr/bin/env python3

"""
.. module:: testUpperLimit
   :synopsis: Test expResultObj.getUpperLimitFor with various inputs.

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.base.physicsUnits import GeV, pb
from databaseLoader import database
from unitTestHelpers import theorySMSFromString as fromString
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.base.model import Model

slhafile = './testFiles/slha/lightEWinos.slha'
model = Model( BSMparticles=BSMList, SMparticles=SMList)
model.updateParticles(inputFile=slhafile,ignorePromptQNumbers=['spin'])


class UpperLimitTest(unittest.TestCase): 

    def testDirectDecay(self):
        expRes=database.getExpResults(analysisIDs = [ "ATLAS-SUSY-2013-05" ], 
                    datasetIDs= [ None ] , txnames= [ "T2bb" ] )[0]
        
        sms = fromString("(PV > sb_1(1),sb_1(2)), (sb_1(1) > b,N1), (sb_1(2) > b,N1)",
                         model=model)
        b1 = model.getParticle(label='sb_1')
        n1 = model.getParticle(label='N1')
        b1.mass = 400.0*GeV
        n1.mass = 100.0*GeV

        txname = expRes.getTxNames()[0]
        smsMatch = txname.hasSMSas(sms)

        ul = expRes.getUpperLimitFor(txname= "T2bb",sms = smsMatch).asNumber(pb)
        self.assertAlmostEqual(ul, 0.0608693)

    def testOutofBounds(self):
        expRes=database.getExpResults(analysisIDs = [ "ATLAS-SUSY-2013-05" ],
                datasetIDs= [ None ] , txnames= [ "T6bbWW" ] )[0]

        sms = fromString("(PV > sb_1(1),sb_1(2)), (sb_1(1) > b,C1+(3)), (sb_1(2) > b,C1+(4)), (C1+(3) > W+,N1), (C1+(4) > W+,N1)",
                         model=model)
        b1 = model.getParticle(label='sb_1')
        c1 = model.getParticle(label='C1+')
        n1 = model.getParticle(label='N1')
        b1.mass = 400.0*GeV
        c1.mass = 250.0*GeV
        n1.mass = 100.0*GeV

        txname = expRes.getTxNames()[0]
        smsMatch = txname.hasSMSas(sms)

        ul = expRes.getUpperLimitFor (txname= "T6bbWW", sms=smsMatch )
        self.assertTrue( ul is None )
 
    def testCascadeDecay(self):

        expRes=database.getExpResults(analysisIDs = [ "ATLAS-SUSY-2013-05" ],
                    datasetIDs= [ None ] , txnames= [ "T6bbWW" ] )[0]

        sms = fromString("(PV > sb_1(1),sb_1(2)), (sb_1(1) > b,C1+(3)), (sb_1(2) > b,C1+(4)), (C1+(3) > W+,N1), (C1+(4) > W+,N1)",
                         model=model)
        b1 = model.getParticle(label='sb_1')
        c1 = model.getParticle(label='C1+')
        n1 = model.getParticle(label='N1')
        b1.mass = 150.0*GeV
        c1.mass = 140.0*GeV
        n1.mass = 135.0*GeV

        txname = expRes.getTxNames()[0]
        smsMatch = txname.hasSMSas(sms)


        ul = expRes.getUpperLimitFor(txname= "T6bbWW",sms=smsMatch).asNumber(pb)
        self.assertAlmostEqual( ul, 324.682 )

if __name__ == "__main__":
    unittest.main()

