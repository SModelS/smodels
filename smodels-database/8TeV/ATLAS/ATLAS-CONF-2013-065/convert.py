#!/usr/bin/env python

"""
.. module:: convert
   :synopsis: used to create info.txt and the <txname>.txt files.

"""
import sys
import os
import argparse

argparser = argparse.ArgumentParser(description =  
'create info.txt, txname.txt, twiki.txt and sms.py')
argparser.add_argument ('-utilsPath', '--utilsPath', 
help = 'path to the package smodels_utils',\
type = str )
argparser.add_argument ('-smodelsPath', '--smodelsPath', 
help = 'path to the package smodels_utils',\
type = str )
args = argparser.parse_args()

if args.utilsPath:
    utilsPath = args.utilsPath
else:
    databaseRoot = '../../../'
    sys.path.append(os.path.abspath(databaseRoot))
    from utilsPath import utilsPath
    utilsPath = databaseRoot + utilsPath
if args.smodelsPath:
    sys.path.append(os.path.abspath(args.smodelsPath))

sys.path.append(os.path.abspath(utilsPath))
from smodels_utils.dataPreparation.inputObjects import MetaInfoInput,DataSetInput
from smodels_utils.dataPreparation.databaseCreation import databaseCreator
from smodels_utils.dataPreparation.massPlaneObjects import x, y, z



#+++++++ global info block ++++++++++++++
info = MetaInfoInput('ATLAS-CONF-2013-065')
info.comment = 'T6bbWWM1300 combines the T6bbWW result from ATLAS-CONF-2013-065 and from ATLAS-CONF-2013-048'
info.sqrts = '8.0'
info.private = False
info.lumi = '20.3'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-065/ http://cds.cern.ch/record/1562840'
info.supersededBy =  "ATLAS-SUSY-2013-19"
info.prettyName = '2 leptons + (b)jets + Etmiss'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.checked ="VM"
T2tt.constraint ="[[['t+']],[['t-']]]"
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T2tt = T2tt.addMassPlane(2*[[x, y]])
T2tt.figure = 'Fig. 7a'
T2tt.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-065/fig_07a.png'
T2tt.dataUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-065/fig_07a_PRELIMINARY.data'
T2tt.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T2tt.txt', 'orig/T2ttOF.data'],
                 dataFormats= ['txt', 'txt'])

#+++++++ next txName block ++++++++++++++
T6bbWW = dataset.addTxName('T6bbWW')
T6bbWW.checked ="VM"
T6bbWW.constraint ="[[['b'],['W+']],[['b'],['W-']]]"
T6bbWW.conditionDescription ="None"
T6bbWW.condition ="None"
T6bbWW.source = "ATLAS"
T6bbWW.massConstraint = None
T6bbWWoff = dataset.addTxName('T6bbWWoff')
T6bbWWoff.constraint ="22*([[['b'],['l+','nu']],[['b'],['l-','nu']]])"
T6bbWWoff.conditionDescription="[[['b'],['l+','nu']],[['b'],['l-','nu']]] > 2*[[['b'],['e+','nu']],[['b'],['e-','nu']]]"
T6bbWWoff.condition="Cgtr([[['b'],['l+','nu']],[['b'],['l-','nu']]],2*[[['b'],['e+','nu']],[['b'],['e-','nu']]])"
T6bbWWoff.massConstraint = [['dm >= 0.0', 'dm <= 76.0'], ['dm >= 0.0', 'dm <= 76.0']]
T6bbWWoff.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T6bbWWM1300 = T6bbWW.addMassPlane(2*[[300.0, x, y]])
T6bbWWM1300.figure = 'Fig. 11'
T6bbWWM1300.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-065/fig_11.png'
T6bbWWM1300.dataUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-065/fig_11_PRELIMINARY.data'
T6bbWWM1300.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T6bbWWM1300.txt', 'orig/T6bbWWCombination.data'],
                 dataFormats= ['txt', 'txt'])
T6bbWWoff.addMassPlane(T6bbWWM1300)



databaseCreator.create()
