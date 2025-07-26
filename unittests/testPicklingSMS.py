#!/usr/bin/env python3

"""
.. module:: testPicklingSMS
   :synopsis: Tests SMS pickling.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import sys
sys.path.insert(0,"../")
import unittest
from smodels.base.smodelsLogging import setLogLevel
from smodels.experiment.expSMS import ExpSMS
from smodels.experiment.defaultFinalStates import finalStates

import pickle as serializer
import os

setLogLevel('error')


class PickilingTest(unittest.TestCase):

    def testPickling(self):
        
        smsString = "[[['jet']],[['jet', 'e-'],['nu','e-']]]"
        sms = ExpSMS.from_string(smsString,model=finalStates)
        
        self.assertEqual([str(n) for n in sms.nodes],
                         ['PV', 'anyBSM', 'anyBSM', 'MET', 'jet', 'anyBSM', 'e-', 'jet', 'MET', 'e-', 'nu'])
        edges = [(str(e1),str(e2)) for e1,e2 in sms.edges]
        self.assertEqual(edges,[('PV', 'anyBSM'), ('PV', 'anyBSM'), ('anyBSM', 'MET'), ('anyBSM', 'jet'), ('anyBSM', 'anyBSM'), ('anyBSM', 'e-'), ('anyBSM', 'jet'), ('anyBSM', 'MET'), ('anyBSM', 'e-'), ('anyBSM', 'nu')])

        testFile = "test_element.pcl"
        with open(testFile,"wb") as f:
            serializer.dump(sms,f,protocol=4)

        with open(testFile,"rb") as f:
            newSMS = serializer.load(f)

        self.assertEqual([str(n) for n in newSMS.nodes],
                            ['PV', 'anyBSM', 'anyBSM', 'MET', 'jet', 'anyBSM', 'e-', 'jet', 'MET', 'e-', 'nu'])
        edges = [(str(e1),str(e2)) for e1,e2 in newSMS.edges]
        self.assertEqual(edges,[('PV', 'anyBSM'), ('PV', 'anyBSM'), ('anyBSM', 'MET'), ('anyBSM', 'jet'), ('anyBSM', 'anyBSM'), ('anyBSM', 'e-'), ('anyBSM', 'jet'), ('anyBSM', 'MET'), ('anyBSM', 'e-'), ('anyBSM', 'nu')])

        os.remove(testFile)

if __name__ == "__main__":
    unittest.main()
