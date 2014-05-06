#!/usr/bin/env python

"""
.. module:: testSlhaChecks
   :synopsis: Tests the slha checker

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import unittest
import setPath

class SlhaTest(unittest.TestCase):
    def testFirstFile(self):
        from smodels.SModelS import installDirectory
        from smodels.tools import slhaChecks

        filename = "%sinputFiles/slha/andrePT4.slha" % (installDirectory() )
        st=slhaChecks.SlhaStatus(filename)
        self.assertEquals ( st.status, (1, 'Input file ok') )
        
        filename = "%sinputFiles/slha/nobdecay.slha" % (installDirectory() )
        st=slhaChecks.SlhaStatus(filename)
        self.assertEquals (st.status, (-1, '#ERROR: charged lsp or long-lived NLSP.\n#Warnings:\n#Empty decay block for PIDs, 1000005.\n#XSECTION table missing.\n#Charged NLSP is stable.\n'))

if __name__ == "__main__":
    unittest.main()
