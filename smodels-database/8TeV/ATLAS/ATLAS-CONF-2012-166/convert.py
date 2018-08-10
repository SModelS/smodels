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
info = MetaInfoInput('ATLAS-CONF-2012-166')
info.comment = 'superseding publication contain more Data'
info.sqrts = '8.0'
info.private = False
info.lumi = '13.0'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2012-166/'
info.supersededBy = 'ATLAS-SUSY-2013-15'
info.prettyName = '1 lepton + 4(1 b-)jets + Etmiss'
info.implementedBy = 'MT'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.checked ="AL"
T2tt.constraint ="[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T2tt = T2tt.addMassPlane(2*[[x, y]])
T2tt.figure = 'Fig.(aux) 3'
T2tt.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2012-166//figaux_03.png'
T2tt.dataUrl = 'Not defined'
T2tt.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T2tt.txt', 'orig/T2tt.txt'],
                 dataFormats= ['txt', 'txt'])

#+++++++ next txName block ++++++++++++++
T6bbWW = dataset.addTxName('T6bbWW')
T6bbWW.checked ="AL"
T6bbWW.constraint ="[[['b'],['W']],[['b'],['W']]]"
T6bbWW.conditionDescription ="None"
T6bbWW.condition ="None"
T6bbWW.source = "ATLAS"
T6bbWW.massConstraint = None
T6bbWWoff = dataset.addTxName('T6bbWWoff')
T6bbWWoff.constraint ="2.3*([[['b'],['L','nu']],[['b'],['jet','jet']]])"
T6bbWWoff.conditionDescription ="[[['b','L','nu']],[['b','jet','jet']]] > 2.7* [[['b','ta','nu']],[['b','jet','jet']]], [[['b','L','nu']],[['b','jet','jet']]] > 2.7* [[['b','e','nu']],[['b','jet','jet']]]"
T6bbWWoff.condition ="Cgtr([[['b'],['L','nu']],[['b'],['jet','jet']]],3.*[[['b'],['ta','nu']],[['b'],['jet','jet']]]);Cgtr([[['b'],['L','nu']],[['b'],['jet','jet']]],3.*[[['b'],['e','nu']],[['b'],['jet','jet']]])"
T6bbWWoff.massConstraint = [['dm >= 0.0', 'dm <= 76.0'], ['dm >= 0.0', 'dm <= 76.0']]
T6bbWWoff.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T6bbWWC150 = T6bbWW.addMassPlane(2*[[x, 150.0, y]])
T6bbWWC150.figure = 'Fig.(aux) 4'
T6bbWWC150.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2012-166//figaux_04.png'
T6bbWWC150.dataUrl = 'Not defined'
T6bbWWC150.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T6bbWWC150_excl.txt', 'orig/T6bbWWC150.txt'],
                 dataFormats= ['txt', 'txt'])
T6bbWWoff.addMassPlane(T6bbWWC150)
#+++++++ next mass plane block ++++++++++++++
T6bbWWx200 = T6bbWW.addMassPlane(2*[[x, y*2.0, y]])
T6bbWWx200.figure = 'Fig.(aux) 5'
T6bbWWx200.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2012-166//figaux_05.png'
T6bbWWx200.dataUrl = 'Not defined'
T6bbWWx200.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T6bbWWx2_excl.txt', 'orig/T6bbWWx2.txt'],
                 dataFormats= ['txt', 'txt'])
T6bbWWoff.addMassPlane(T6bbWWx200)



databaseCreator.create()
