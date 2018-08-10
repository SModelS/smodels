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
info = MetaInfoInput('ATLAS-CONF-2013-037')
info.sqrts = '8.0'
info.private = False
info.lumi = '20.7'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-037/'
info.supersededBy = 'ATLAS-SUSY-2013-15'
info.contact = 'Sophio Pataria <Sophio.Pataraia@cern.ch>, Tobias Golling <tobias.golling@cern.ch>'
info.prettyName = '1 lepton + >= 4(1 b-)jets + Etmiss'
# info.supersedes = 'ATLAS-CONF-2012-166'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.constraint ="[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T2tt = T2tt.addMassPlane(2*[[x, y]])
T2tt.figure = 'Fig. 4'
T2tt.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-037/figaux_04.png'
T2tt.dataUrl = 'Not defined'
T2tt.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T2tt.txt', 'orig/T2tt_corr.txt'],
                 dataFormats= ['txt', 'txt'])

#+++++++ next txName block ++++++++++++++
T6bbWW = dataset.addTxName('T6bbWW')
T6bbWW.constraint ="[[['b'],['W']],[['b'],['W']]]"
T6bbWW.conditionDescription ="None"
T6bbWW.condition ="None"
T6bbWW.source = "ATLAS"
T6bbWW.massConstraint = None
T6bbWWoff = dataset.addTxName('T6bbWWoff')
T6bbWWoff.constraint ="3.5*([[['b'],['L','nu']],[['b'],['jet','jet']]])"
T6bbWWoff.conditionDescription ="[[['b','L','nu']],[['b','jet','jet']]] > 2.7* [[['b','ta','nu']],[['b','jet','jet']]], [[['b','L','nu']],[['b','jet','jet']]] > 2.7* [[['b','e','nu']],[['b','jet','jet']]]"
T6bbWWoff.condition ="Cgtr([[['b'],['L','nu']],[['b'],['jet','jet']]],3.*[[['b'],['ta','nu']],[['b'],['jet','jet']]]);Cgtr([[['b'],['L','nu']],[['b'],['jet','jet']]],3.*[[['b'],['e','nu']],[['b'],['jet','jet']]])"
T6bbWWoff.massConstraint = [['dm >= 0.0', 'dm <= 76.0'], ['dm >= 0.0', 'dm <= 76.0']]
T6bbWWoff.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T6bbWWC150 = T6bbWW.addMassPlane(2*[[x, 150.0, y]])
T6bbWWC150.figure = 'Fig.(aux) 5'
T6bbWWC150.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-037/figaux_05.png'
T6bbWWC150.dataUrl = 'Not defined'
T6bbWWC150.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T6bbWWC150.txt', 'orig/T6bbWWC150_corr.txt'],
                 dataFormats= ['txt', 'txt'])
T6bbWWoff.addMassPlane(T6bbWWC150)
#+++++++ next mass plane block ++++++++++++++
T6bbWWx200 = T6bbWW.addMassPlane(2*[[x, y*2.0, y]])
T6bbWWx200.figure = 'Fig.(aux) 6'
T6bbWWx200.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-037/figaux_06.png'
T6bbWWx200.dataUrl = 'Not defined'
T6bbWWx200.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T6bbWWx200.txt', 'orig/T6bbWWx200_corr.txt'],
                 dataFormats= ['txt', 'txt'])
T6bbWWoff.addMassPlane(T6bbWWx200)



databaseCreator.create()
