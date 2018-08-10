#!/usr/bin/env python

"""
.. module:: convert
   :synopsis: used to create info.txt and the <txname>.txt files.

"""
import sys
import os
import argparse
import types

argparser = argparse.ArgumentParser(description =  
'create info.txt, txname.txt, twiki.txt and sms.py')
argparser.add_argument ('-utilsPath', '--utilsPath', 
help = 'path to the package smodels_utils',\
type = str)
argparser.add_argument ('-smodelsPath', '--smodelsPath', 
help = 'path to the package smodels_utils',\
type = str)
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
info = MetaInfoInput('ATLAS-CONF-2013-001')
info.sqrts = '8.0'
info.private = False
info.lumi = '12.8'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-001/'
info.supersededBy = 'ATLAS-SUSY-2013-05'
info.prettyName = '0 leptons + 2 b-jets + Etmiss'
info.implementedBy = 'MT'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T6bbWW = dataset.addTxName('T6bbWW')
T6bbWW.constraint = "[[['b'],['W']],[['b'],['W']]]"
T6bbWW.condition = "None"
T6bbWW.massConstraint = [['dm >= 0.0', 'dm >= 76.0'], ['dm >= 0.0', 'dm >= 76.0']]
T6bbWW.source = 'ATLAS'
T6bbWWoff = dataset.addTxName('T6bbWWoff')
T6bbWWoff.checked ="AL"
T6bbWWoff.constraint ="[[['b'],['L','nu']],[['b'],['L','nu']]] + [[['b'],['L','nu']],[['b'],['jet','jet']]] + [[['b'],['jet','jet']],[['b'],['jet','jet']]]"
T6bbWWoff.conditionDescription ="None"
T6bbWWoff.condition ="None"
T6bbWWoff.massConstraint = [['dm >= 0.0', 'dm <= 76.0'], ['dm >= 0.0', 'dm <= 76.0']]
T6bbWWoff.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T6bbWWD020 = T6bbWW.addMassPlane(2*[[x, y+20.0, y]])
T6bbWWD020.figure ='Fig.(aux) 7'
T6bbWWD020.figureUrl ='https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-001/figaux_07.png'
T6bbWWD020.dataUrl = 'Not defined'
T6bbWWD020.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T6bbWW20GeV_excl.txt', 'orig/T6bbWWoff.txt'],
                 dataFormats= ['txt', 'txt'])
T6bbWWoff.addMassPlane(T6bbWWD020)
#+++++++ next mass plane block ++++++++++++++
T6bbWWD005 = T6bbWW.addMassPlane(2*[[x, y+5.0, y]])
T6bbWWD005.figure ='Fig.(aux) 6'
T6bbWWD005.figureUrl ='https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-001/figaux_06.png'
T6bbWWD005.dataUrl = 'Not defined'
T6bbWWD005.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T6bbWW5GeV_excl.txt', 'orig/T2bb.txt'],
                 dataFormats= ['txt', 'txt'])
T6bbWWoff.addMassPlane(T6bbWWD005)



databaseCreator.create()
