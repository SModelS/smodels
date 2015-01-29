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

class InterpolationTest(unittest.TestCase):
    def testInterpolation(self):
        from smodels.experiment import infoObjects, dataObjects
        from smodels.tools.physicsUnits import GeV, TeV
        path="/home/walten/Downloads/extended-database/"
        expid="8TeV/CMS/CMS-SUS-12-028/"
        info = infoObjects.InfoFile(os.path.join(path,expid,"info.txt"))
        data = dataObjects.DataFile(os.path.join(path,expid,"sms.py"),info)
        for i in data.dataList:
            if i.txname != "T1":
                continue
            print i.analysisID, i.txname
            ## print i.data
