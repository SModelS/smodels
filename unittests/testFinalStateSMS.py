#!/usr/bin/env python3

"""
.. module:: testElementClass
   :synopsis: Tests the theory.element.Element class

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""
import sys
sys.path.insert(0,"../")
import unittest
from smodels.share.models.mssm import BSMList
from smodels.share.models.SMparticles import SMList
from smodels.base.model import Model
from collections import OrderedDict
from smodels.base.physicsUnits import fb, GeV, TeV
from smodels.share.models import SMparticles
from smodels.experiment.defaultFinalStates import finalStates
from unitTestHelpers import theorySMSFromString as fromString
from smodels.tools.coverage import FinalStateSMS
from smodels.experiment.defaultFinalStates import lList, tList, MET




slhafile = './testFiles/slha/lightEWinos.slha'

model = Model( BSMparticles=BSMList, SMparticles=SMList)
model.updateParticles(inputFile=slhafile,ignorePromptQNumbers=['spin'])

class FinalStateCompTest(unittest.TestCase):

    def testFinalState(self):


        stringEl = "(PV > gluino(1),gluino(2)), (gluino(1) > b,b,N2(3)), (gluino(2) > b,t+,C1-(4)), (N2(3) > N1,b,b), (C1-(4) > N1~,q,q)"
        sms = fromString(stringEl, model=model)

        # Compress SMS:
        fsSMS = FinalStateSMS(sms,smFinalStates=[lList,tList],bsmFinalStates=[MET],
                              missingX=10*fb)
        compStr = 'PV > (b,t,q,q,MET), (b,b,b,b,MET)'
        self.assertEqual(str(fsSMS),compStr)
        self.assertEqual(fsSMS.missingX,10*fb)


        stringEl = "(PV > gluino(1),gluino(2)), (gluino(1) > mu-,nu,b,N2(3)), (gluino(2) > b,t+,C1-(4)), (N2(3) > N1,b,b), (C1-(4) > N1~,e-,e+)"
        sms = fromString(stringEl, model=model)

        # Compress SMS:
        fsSMS = FinalStateSMS(sms,smFinalStates=[lList,tList],bsmFinalStates=[MET],
                              missingX=1.5*fb)
        compStr = 'PV > (b,t,l,l,MET), (b,b,b,nu,l,MET)'
        self.assertEqual(str(fsSMS),compStr)
        self.assertEqual(fsSMS.missingX,1.5*fb)


if __name__ == "__main__":
    unittest.main()
