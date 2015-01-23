#!/usr/bin/env python

"""
.. module:: Example
   :synopsis: Basic main file example for using SModelS.

"""

#Import basic functions (this file must be run under the installation folder)
import sys
from smodels.experiment.infoObjects import InfoFile
from smodels.experiment.dataObjects import DataFile
from smodels.experiment.databaseObjects import DataBase
from smodels.tools.physicsUnits import GeV, fb, TeV, pb
from smodels.experiment.databaseBrowser import Browser
import numpy as np

database = DataBase("/home/lessa/smodels-DBtest/")
browser = Browser(database)
#print browser.getAttributes()
#print browser.getAttributes()
nl = 0
print browser
browser.loadExpResultsWith({'txname': ['T2','T1','T2bb'], 'id' : ['ATLAS-SUSY-2013-05']})
print browser
print browser.getValuesFor("implimented_by")
print browser.getValuesFor("constraint")

#sys.exit()
database = DataBase("/home/lessa/smodels-DBtest/")
print database._getDatabaseVersion
res = database._getExpResults(txnames=["T2bb"])
for r in res: print r.info.getInfo('id'),r.info.getTxNames()
for r in res: print r.info.getInfo('constraint')
for r in res: print r.info.getInfo('publication')
sys.exit()

listOfana = database.getAnalyses(txnames=['T2bb'])


for ana in listOfana:
  print ana.printout()
#  print ana.getUpperLimitFor([[200.*GeV,50.*GeV],[200.*GeV,50.*GeV]])
