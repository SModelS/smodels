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
info = MetaInfoInput('ATLAS-SUSY-2013-15')
info.sqrts = '8.0'
info.private = False
info.lumi = '20.3'
info.publication =  'link.springer.com/article/10.1007/JHEP11(2014)118'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-15/'
info.arxiv = 'http://arxiv.org/abs/1407.0583'
info.contact =  "ATLAS collaboration"
info.prettyName = '1 lepton + 4 (1 b-)jets + Etmiss'
info.supersedes = 'ATLAS-CONF-2013-037; ATLAS-CONF-2012-166'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T2bbWW = dataset.addTxName('T2bbWW')
T2bbWW.constraint ="[[['b','W']],[['b','W']]]"
T2bbWW.conditionDescription ="None"
T2bbWW.condition ="None"
T2bbWW.source = "ATLAS"
T2bbWW.massConstraint = None
T2bbWWoff = dataset.addTxName('T2bbWWoff')
T2bbWWoff.constraint ="2.3*([[['b','L','nu']],[['b','jet','jet']]])"
T2bbWWoff.conditionDescription ="[[['b','L','nu']],[['b','jet','jet']]] > 2.7* [[['b','ta','nu']],[['b','jet','jet']]], [[['b','L','nu']],[['b','jet','jet']]] > 2.7* [[['b','e','nu']],[['b','jet','jet']]]"
T2bbWWoff.condition ="Cgtr([[['b','L','nu']],[['b','jet','jet']]],3.*[[['b','ta','nu']],[['b','jet','jet']]]);Cgtr([[['b','L','nu']],[['b','jet','jet']]],3.*[[['b','e','nu']],[['b','jet','jet']]])"
T2bbWWoff.massConstraint = [['dm <= 76.0'], ['dm <= 76.0']]
T2bbWWoff.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T2bbWW = T2bbWW.addMassPlane(2*[[x, y]])
T2bbWW.figure = 'Fig.(aux) 14'
T2bbWW.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-15/figaux_014.png'
T2bbWW.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1304456/d4'
T2bbWW.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T2bbWW.txt', 'orig/limit_T2bbWW.txt'],
                 dataFormats= ['txt', 'txt'])
T2bbWWoff.addMassPlane(T2bbWW)

#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.constraint ="[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T2tt = T2tt.addMassPlane(2*[[x, y]])
T2tt.figure = 'Fig.(aux) 13'
T2tt.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-15/figaux_013.png'
T2tt.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1304456/d1'
T2tt.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T2tt.txt', 'orig/limit_T2tt.txt'],
                 dataFormats= ['txt', 'txt'])

#+++++++ next txName block ++++++++++++++
T6bbWW = dataset.addTxName('T6bbWW')
T6bbWW.constraint ="[[['b'],['W']],[['b'],['W']]]"
T6bbWW.conditionDescription = "None"
T6bbWW.condition = "None"
T6bbWW.source = "ATLAS"
T6bbWW.massConstraint = None
T6bbWWoff = dataset.addTxName('T6bbWWoff')
T6bbWWoff.constraint ="2.3*([[['b'],['L','nu']],[['b'],['jet','jet']]])"
T6bbWWoff.conditionDescription ="[[['b'],['L','nu']],[['b'],['jet','jet']]] > 2.7* [[['b'],['ta','nu']],[['b'],['jet','jet']]],[[['b'],['L','nu']],[['b'],['jet','jet']]] > 2.7* [[['b'],['e','nu']],[['b'],['jet','jet']]]"
T6bbWWoff.condition ="Cgtr([[['b'],['L','nu']],[['b'],['jet','jet']]],3.*[[['b'],['ta','nu']],[['b'],['jet','jet']]]);Cgtr([[['b'],['L','nu']],[['b'],['jet','jet']]],3.*[[['b'],['e','nu']],[['b'],['jet','jet']]])"
T6bbWWoff.massConstraint = [['dm >= 0.0', 'dm <= 76.0'], ['dm >= 0.0', 'dm <= 76.0']]
T6bbWWoff.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T6bbWWoffD10 = T6bbWW.addMassPlane(2*[[x, x-10., y]])
T6bbWWoffD10.figure = "fig(aux) 21"
T6bbWWoffD10.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-15/figaux_021.png"
T6bbWWoffD10.dataUrl = 'Not defined'
T6bbWWoffD10.setSources(dataLabels= ['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_exp_T6bbWWoffD010.txt', 'orig/exclusion_T6bbWWoffD010.txt', 'orig/limit_T6bbWWoffD010.txt'],
                 dataFormats= ['txt', 'txt', 'txt'])
T6bbWWoff.addMassPlane(T6bbWWoffD10)
#+++++++ next mass plane block ++++++++++++++
T6bbWWoffD5 = T6bbWWoff.addMassPlane(2*[[x, y+5., y]]) #This is a purely off-shell plane, do not add it to the on-shell txname
T6bbWWoffD5.figure = "fig(aux) 19"
T6bbWWoffD5.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-15/figaux_019.png"
T6bbWWoffD5.dataUrl = 'Not defined'
T6bbWWoffD5.setSources(dataLabels= ['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_exp_T6bbWWoffD005.txt', 'orig/exclusion_T6bbWWoffD005.txt', 'orig/limit_T6bbWWoffD005.txt'],
                 dataFormats= ['txt', 'txt', 'txt'])
#+++++++ next mass plane block ++++++++++++++
T6bbWWoffD20 = T6bbWWoff.addMassPlane(2*[[x, y+20., y]]) #This is a purely off-shell plane, do not add it to the on-shell txname
T6bbWWoffD20.figure = "fig(aux) 20"
T6bbWWoffD20.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-15/figaux_020.png"
T6bbWWoffD20.dataUrl = 'Not defined'
T6bbWWoffD20.setSources(dataLabels= ['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_exp_T6bbWWoffD020.txt', 'orig/exclusion_T6bbWWoffD020.txt', 'orig/limit_T6bbWWoffD020.txt'],
                 dataFormats= ['txt', 'txt', 'txt'])
#+++++++ next mass plane block ++++++++++++++
T6bbWWoffM1300 = T6bbWW.addMassPlane(2*[[300, x, y]])
T6bbWWoffM1300.figure = "fig(aux) 22"
T6bbWWoffM1300.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-15/figaux_022.png"
T6bbWWoffM1300.dataUrl = 'Not defined'
T6bbWWoffM1300.setSources(dataLabels= ['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_exp_T6bbWWoffM1300.txt', 'orig/exclusion_T6bbWWoffM1300.txt', 'orig/limit_T6bbWWoffM1300.txt'],
                 dataFormats= ['txt', 'txt', 'txt'])
T6bbWWoff.addMassPlane(T6bbWWoffM1300)
#+++++++ next mass plane block ++++++++++++++
T6bbWWoffx200 = T6bbWW.addMassPlane(2*[[x, y*2.0, y]])
T6bbWWoffx200.figure = "fig(aux) 16"
T6bbWWoffx200.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-15/figaux_016.png"
T6bbWWoffx200.dataUrl = 'Not defined'
T6bbWWoffx200.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T6bbWWoffx200.txt', 'orig/limit_T6bbWWoffx200.txt'],
                 dataFormats= ['txt', 'txt'])
T6bbWWoff.addMassPlane(T6bbWWoffx200)
#+++++++ next mass plane block ++++++++++++++
T6bbWWoffC106 = T6bbWW.addMassPlane(2*[[x, 106.0, y]])
T6bbWWoffC106.figure = "fig(aux) 18"
T6bbWWoffC106.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-15/figaux_018.png"
T6bbWWoffC106.dataUrl = 'Not defined'
T6bbWWoffC106.setSources(dataLabels= ['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_exp_T6bbWWoffC106.txt', 'orig/exclusion_T6bbWWoffC106.txt', 'orig/limit_T6bbWWoffC106.txt'],
                 dataFormats= ['txt', 'txt', 'txt'])
T6bbWWoff.addMassPlane(T6bbWWoffC106)
#+++++++ next mass plane block ++++++++++++++
T6bbWWoffC150 = T6bbWW.addMassPlane(2*[[x, 150.0, y]])
T6bbWWoffC150.figure = "fig(aux) 17"
T6bbWWoffC150.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-15/figaux_017.png"
T6bbWWoffC150.dataUrl = 'Not defined'
T6bbWWoffC150.setSources(dataLabels= ['expExclusion', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_exp_T6bbWWoffC150.txt', 'orig/exclusion_T6bbWWoffC150.txt', 'orig/limit_T6bbWWoffC150.txt'],
                 dataFormats= ['txt', 'txt', 'txt'])
T6bbWWoff.addMassPlane(T6bbWWoffC150)



databaseCreator.create()
