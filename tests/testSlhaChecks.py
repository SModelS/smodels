#!/usr/bin/env python

"""
.. module:: testSlhaChecks
   :synopsis: Tests the slha checker

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>

"""
import unittest
import setPath
from smodels.installation import installDirectory
from smodels.tools import slhaChecks

class SlhaTest(unittest.TestCase):
    def testGoodFile(self):

        filename = "%sinputFiles/slha/andrePT4.slha" % (installDirectory() )
        st=slhaChecks.SlhaStatus(filename)
        self.assertEquals ( st.status, (1, 'Input file ok') )
        
    def testBadFile(self):
        filename = "%sinputFiles/slha/nobdecay.slha" % (installDirectory() )
        st=slhaChecks.SlhaStatus(filename)
        self.assertEquals (st.status, (-1, '#ERROR: special signatures in this point.\n#Warnings:\n#Empty decay block for PIDs, 1000005.\n#XSECTION table missing, will be computed by SModelS.\n#Charged NLSP is stable.\n'))

if __name__ == "__main__":
    unittest.main()
