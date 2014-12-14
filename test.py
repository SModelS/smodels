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

database = DataBase("/home/lessa/smodels-HSCP/smodels-database/")
browser = Browser(database)
#print browser.getAttributes()
#print browser.getAttributes()
nl = 0
print browser
browser.loadExpResultsWith({'txname': ['T2','T1'], 'id' : ['ATLAS-CONF-2013-047','CMS-SUS-12-028','CMS-SUS-12-018']})
print browser

sys.exit()
database = DataBase("/home/lessa/smodels-HSCP/smodels-database/")
print database._getDatabaseVersion
res = database._getExpResults(txnames=["TChiWZ"])
for r in res: print r.info.getInfo('id'),r.info.getTxNames()

listOfana = database.getAnalyses()


for ana in listOfana:
  print ana.printout()
#  print ana.getUpperLimitFor([[200.*GeV,50.*GeV],[200.*GeV,50.*GeV]])
