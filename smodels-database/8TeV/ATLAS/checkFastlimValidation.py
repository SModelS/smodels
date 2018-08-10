#!/usr/bin/env python


import sys
sys.path.append('/home/lessa/smodels/')
sys.path.append('/home/lessa/smodels-utils/')
from smodels.experiment.databaseObj import Database
from smodels.tools.databaseBrowser import Browser


database = Database('/home/lessa/smodels-database/')

allResults = [exp  for exp in database.getExpResults(useSuperseded = True, useNonValidated = True)]
fastlimResults = {}
for exp in allResults:
    if hasattr(exp.globalInfo,'contact') and exp.globalInfo.contact == 'fastlim':
        fastlimResults[exp.globalInfo.id] = exp
print len(fastlimResults)

from validation.plottingFuncs import getExclusionCurvesFor
canBeValidated = {}
for exp in allResults:
    if exp.globalInfo.id in fastlimResults:
        canBeValidated[exp.globalInfo.id] = []
        for tx in fastlimResults[exp.globalInfo.id].getTxNames():
            tgraphs = getExclusionCurvesFor(fastlimResults[exp.globalInfo.id],tx.txName)
            if tgraphs:
                canBeValidated[exp.globalInfo.id].append({'txname' : tx, 'graphs': tgraphs})


for  exp in canBeValidated:
    if  not canBeValidated[exp]: continue
    print exp
    txs = set([tx['txname'].txName for tx in canBeValidated[exp]])
    print list(txs)

doneValidationPlot = {'ATLAS-CONF-2013-048' : '17/05/2016', 'ATLAS-CONF-2013-053' : '17/05/2016', 
                      'ATLAS-CONF-2013-061' : '22/05/2016'}

