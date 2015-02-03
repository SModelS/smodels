#!/usr/bin/env python

"""
.. module:: Example
   :synopsis: Basic main file example for using SModelS.

"""

#Import basic functions (this file must be run under the installation folder)
import sys
from smodels.experiment.infoObject import Info
from smodels.experiment.databaseObjects import DataBase
from smodels.tools.physicsUnits import GeV, fb, TeV, pb
from smodels.experiment.databaseBrowser import Browser
import numpy as np

# database = DataBase("/home/lessa/smodels-database/")
# browser = Browser(database)
# print browser.getAttributes()
#print browser.getAttributes()
# nl = 0
# print browser
# browser.loadExpResultsWith({'txname': ['T2','T1','T2bb'], 'id' : ['ATLAS-SUSY-2013-05']})
# print browser
# print browser.getValuesFor("lumi")
# print browser.getValuesFor("axes")
# print database
database = DataBase("/home/lessa/smodels-database/")
print database
listOfExpRes = database.getExpResults()
for expRes in listOfExpRes:
    for txname in expRes.txnames:
        print txname.axes
