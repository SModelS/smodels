#!/usr/bin/env python

"""
.. module:: testInterpolation
   :synopsis: Tests the retrieval of the upper limits, including the PCA.
              Will replace the upper limit test

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""
import unittest
import os
import numpy as np
from smodels.tools.physicsUnits import GeV, TeV

class InterpolationTest(unittest.TestCase):
    def mestInterpolation(self):
        from smodels.experiment import infoObjects, dataObjects
        path="/home/walten/git/smodels-database/" ## new style!
        expid="8TeV/ATLAS/ATLAS-SUSY-2013-05/"
        info = infoObjects.InfoFile(os.path.join(path,expid,"info.txt"))
        data = dataObjects.DataFile(os.path.join(path,expid,"sms.py"),info)
        d=data.getData ( "T2bb" )
        print d.analysisID, d.txname
        ## print i.data
        massarray=[[ 400.*GeV, 100.*GeV ], [ 400.*GeV, 100.*GeV ] ]
        ## print np.array(massarray).shape
        # d.getULFor ( massarray )
    
    def testWithBrowser(self):
        from smodels.experiment import databaseBrowser
        path="/home/walten/git/smodels-database/" ## new style!
        browser = databaseBrowser.Browser ( path )
        txname,id="T2bb","ATLAS-SUSY-2013-05"
        massarray=[[ 400.*GeV, 100.*GeV ], [ 400.*GeV, 100.*GeV ] ]
        browser.getULFor ( id, txname, massarray )
        ## print browser.getAttributes()
        #browser.loadExpResultsWith ( { "txname":txname, "id": id } )
        #for expres in browser:
        #    if expres.info.globalInfo.id != id:
        #        continue
        #    print "ExpRes",expres,expres.info.globalInfo.id
        #    print expres.data.getData ( txname )
        #    print expres.data.getData ( txname ).getULFor ( massarray )

if __name__ == "__main__":
    unittest.main()
