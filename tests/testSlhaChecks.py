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


if __name__ == "__main__":
    unittest.main()
