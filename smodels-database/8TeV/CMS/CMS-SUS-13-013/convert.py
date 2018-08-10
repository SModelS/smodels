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
info = MetaInfoInput('CMS-SUS-13-013')
info.url ='https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13013'
info.sqrts = 8
info.lumi = 19.5
info.prettyName = '2 SS leptons + (b-)jets + Etmiss'
info.private = False
info.arxiv = 'http://arxiv.org/pdf/1311.6736v2.pdf'
info.contact =''
info.publication ='http://link.springer.com/article/10.1007%2FJHEP01%282014%29163'
info.comment =''
info.supersedes ='CMS-SUS-12-017'
info.supersededBy =''


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked =''
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = 'None'
T1tttt.condition = 'None'
T1tttt.source = "CMS"
T1tttt.massConstraint = None
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.constraint ="[[['b','b','W','W']],[['b','b','W','W']]]"
T1ttttoff.conditionDescription ='None'
T1ttttoff.condition ="None"
T1ttttoff.massConstraint = [['dm <= 338.0'], ['dm <= 338.0']]
T1ttttoff.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane(2*[[x, y]])
T1tttt_1.figure ='Figure 5'
T1tttt_1.figureUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/Figure5_A1.pdf'
T1tttt_1.dataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/ModelA1.root'
T1tttt_1.histoDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/Figure5_A1.pdf'
T1tttt_1.dataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/ModelA1.root'
T1tttt_1.exclusionDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/ModelA1.root'
T1tttt_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'expectedUpperLimits', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/ModelA1.root', 'orig/ModelA1.root', 'orig/ModelA1.root', 'orig/ModelA1.root', 'orig/ModelA1.root', 'orig/ModelA1.root', 'orig/ModelA1.root', 'orig/ModelA1.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['graph_smoothed_Exp', 'graph_smoothed_ExpM', 'graph_smoothed_ExpP', 'expectedUL', 'graph_smoothed_Obs', 'graph_smoothed_ObsM', 'graph_smoothed_ObsP', 'observedUL'],units= [None, None, None, 'fb', None, None, None, 'fb'])
T1ttttoff.addMassPlane(T1tttt_1)

#+++++++ next txName block ++++++++++++++
T6ttWW = dataset.addTxName('T6ttWW')
T6ttWW.checked =''
T6ttWW.constraint ="[[['t'],['W']],[['t'],['W']]]"
T6ttWW.conditionDescription = "None"
T6ttWW.condition = "None"
T6ttWW.conditionDescription = None
T6ttWW.condition =None
T6ttWW.source = "CMS"
T6ttWW.massConstraint = None
T6ttWWoff = dataset.addTxName('T6ttWWoff')
T6ttWWoff.checked =''
T6ttWWoff.constraint ="2.3*([[['t'],['L','nu']],[['t'],['jet','jet']]])"
T6ttWWoff.conditionDescription ="[[['t'],['L','nu']],[['t'],['jet','jet']]] > 2.7* [[['t'],['ta','nu']],[['t'],['jet','jet']]],[[['t'],['L','nu']],[['t'],['jet','jet']]] > 2.7* [[['t'],['e','nu']],[['t'],['jet','jet']]]"
T6ttWWoff.condition ="Cgtr([[['t'],['L','nu']],[['t'],['jet','jet']]],3.*[[['t'],['ta','nu']],[['t'],['jet','jet']]]);Cgtr([[['t'],['L','nu']],[['t'],['jet','jet']]],3.*[[['t'],['e','nu']],[['t'],['jet','jet']]])"
T6ttWWoff.massConstraint = [['dm >= 0.0', 'dm <= 76.0'], ['dm >= 0.0', 'dm <= 76.0']]
T6ttWWoff.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T6ttWW_1 = T6ttWW.addMassPlane(2*[[x, y, 50]])
T6ttWW_1.figure = 'Fig.6'
T6ttWW_1.figureUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/Figure6_B1a.pdf'
T6ttWW_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13013'
T6ttWW_1.histoDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/ModelB1.root'
T6ttWW_1.dataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/ModelB1.root'
T6ttWW_1.exclusionDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/ModelB1.root'
T6ttWW_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'expectedUpperLimits', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/ModelB1.root', 'orig/ModelB1.root', 'orig/ModelB1.root', 'orig/ModelB1.root', 'orig/ModelB1.root', 'orig/ModelB1.root', 'orig/ModelB1.root', 'orig/ModelB1.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['graph_smoothed_Exp', 'graph_smoothed_ExpM', 'graph_smoothed_ExpP', 'expectedUL', 'graph_smoothed_Obs', 'graph_smoothed_ObsM', 'graph_smoothed_ObsP', 'observedUL'],units= [None, None, None, 'fb', None, None, None, 'fb'])
T6ttWWoff.addMassPlane(T6ttWW_1)
#+++++++ next mass plane block ++++++++++++++
T6ttWW_2 = T6ttWW.addMassPlane(2*[[x, y/(0.5), y]])
T6ttWW_2.dataUrl = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13013'
T6ttWW_2.histoDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/ModelB1_x0p5.root'
T6ttWW_2.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/ModelB1_x0p5.root'
T6ttWW_2.exclusionDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/ModelB1_x0p5.root'
T6ttWW_2.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'expectedUpperLimits', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/ModelB1_x0p5.root', 'orig/ModelB1_x0p5.root', 'orig/ModelB1_x0p5.root', 'orig/ModelB1_x0p5.root', 'orig/ModelB1_x0p5.root', 'orig/ModelB1_x0p5.root', 'orig/ModelB1_x0p5.root', 'orig/ModelB1_x0p5.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['graph_smoothed_Exp', 'graph_smoothed_ExpM', 'graph_smoothed_ExpP', 'expectedUL', 'graph_smoothed_Obs', 'graph_smoothed_ObsM', 'graph_smoothed_ObsP', 'observedUL'],units= [None, None, None, 'fb', None, None, None, 'fb'])
T6ttWWoff.addMassPlane(T6ttWW_2)
#+++++++ next mass plane block ++++++++++++++
T6ttWW_3 = T6ttWW.addMassPlane(2*[[x, y/(0.8), y]])
T6ttWW_3.dataUrl = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13013'
T6ttWW_3.histoDataUrl =  'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/ModelB1_x0p8.root'
T6ttWW_3.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/ModelB1_x0p8.root'
T6ttWW_3.exclusionDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/ModelB1_x0p8.root'
T6ttWW_3.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'expectedUpperLimits', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/ModelB1_x0p8.root', 'orig/ModelB1_x0p8.root', 'orig/ModelB1_x0p8.root', 'orig/ModelB1_x0p8.root', 'orig/ModelB1_x0p8.root', 'orig/ModelB1_x0p8.root', 'orig/ModelB1_x0p8.root', 'orig/ModelB1_x0p8.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['graph_smoothed_Exp', 'graph_smoothed_ExpM', 'graph_smoothed_ExpP', 'expectedUL', 'graph_smoothed_Obs', 'graph_smoothed_ObsM', 'graph_smoothed_ObsP', 'observedUL'],units= [None, None, None, 'fb', None, None, None, 'fb'])
T6ttWWoff.addMassPlane(T6ttWW_3)



databaseCreator.create()
