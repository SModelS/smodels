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
info.publication ='http://link.springer.com/article/10.1007%2FJHEP01%282014%29163'
info.supersedes ='CMS-SUS-12-017'




#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR27_HighPt")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR27_HighPt", observedN = 0, expectedBG = 1.22 , bgError = 0.62, upperLimit = '1.539E-01*fb', expectedUpperLimit = '2.128E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.constraint = "[[['b','b','W','W']],[['b','b','W','W']]]"
T1ttttoff.conditionDescription ="None"
T1ttttoff.condition = "None"
T1ttttoff.massConstraint = [['dm <= 338.']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure ='Figure 5'
T1tttt_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/Figure5_A1.pdf"
T1tttt_1.addSource('obsExclusion','orig/ModelA1.root', 'root', objectName = 'graph_smoothed_Obs')
T1tttt_1.addSource('obsExclusionM1','orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ObsM')
T1tttt_1.addSource('obsExclusionP1','orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ObsP')
T1tttt_1.addSource('expExclusion', 'orig/ModelA1.root', 'root', objectName = 'graph_smoothed_Exp')
T1tttt_1.addSource('expExclusionM1', 'orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ExpM')
T1tttt_1.addSource('expExclusionP1', 'orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ExpP')
#----limit source----
T1tttt_1.addSource('efficiencyMap', 'orig/SR27_HighPt.txt', 'txt')
T1tttt_1.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/ModelA1.root'
T1ttttoff.addMassPlane(T1tttt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR25_HighPt")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR25_HighPt", observedN = 4, expectedBG = 2.86 , bgError = 1.14, upperLimit = '3.578E-01*fb', expectedUpperLimit = '2.478E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.constraint = "[[['b','b','W','W']],[['b','b','W','W']]]"
T1ttttoff.conditionDescription ="None"
T1ttttoff.condition = "None"
T1ttttoff.massConstraint = [['dm <= 338.']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure ='Figure 5'
T1tttt_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/Figure5_A1.pdf"
T1tttt_1.addSource('obsExclusion','orig/ModelA1.root', 'root', objectName = 'graph_smoothed_Obs')
T1tttt_1.addSource('obsExclusionM1','orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ObsM')
T1tttt_1.addSource('obsExclusionP1','orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ObsP')
T1tttt_1.addSource('expExclusion', 'orig/ModelA1.root', 'root', objectName = 'graph_smoothed_Exp')
T1tttt_1.addSource('expExclusionM1', 'orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ExpM')
T1tttt_1.addSource('expExclusionP1', 'orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ExpP')
#----limit source----
T1tttt_1.addSource('efficiencyMap', 'orig/SR25_HighPt.txt', 'txt')
T1tttt_1.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/ModelA1.root'
T1ttttoff.addMassPlane(T1tttt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR26_HighPt")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR26_HighPt", observedN = 1, expectedBG = 0.81 , bgError = 0.54, upperLimit = '2.188E-01*fb', expectedUpperLimit = '1.540E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.constraint = "[[['b','b','W','W']],[['b','b','W','W']]]"
T1ttttoff.conditionDescription ="None"
T1ttttoff.condition = "None"
T1ttttoff.massConstraint = [['dm <= 338.']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure ='Figure 5'
T1tttt_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/Figure5_A1.pdf"
T1tttt_1.addSource('obsExclusion','orig/ModelA1.root', 'root', objectName = 'graph_smoothed_Obs')
T1tttt_1.addSource('obsExclusionM1','orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ObsM')
T1tttt_1.addSource('obsExclusionP1','orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ObsP')
T1tttt_1.addSource('expExclusion', 'orig/ModelA1.root', 'root', objectName = 'graph_smoothed_Exp')
T1tttt_1.addSource('expExclusionM1', 'orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ExpM')
T1tttt_1.addSource('expExclusionP1', 'orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ExpP')
#----limit source----
T1tttt_1.addSource('efficiencyMap', 'orig/SR26_HighPt.txt', 'txt')
T1tttt_1.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/ModelA1.root'
T1ttttoff.addMassPlane(T1tttt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR23_HighPt")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR23_HighPt", observedN = 3, expectedBG = 3.78 , bgError = 1.44, upperLimit = '2.865E-01*fb', expectedUpperLimit = '2.856E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.constraint = "[[['b','b','W','W']],[['b','b','W','W']]]"
T1ttttoff.conditionDescription ="None"
T1ttttoff.condition = "None"
T1ttttoff.massConstraint = [['dm <= 338.']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure ='Figure 5'
T1tttt_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/Figure5_A1.pdf"
T1tttt_1.addSource('obsExclusion','orig/ModelA1.root', 'root', objectName = 'graph_smoothed_Obs')
T1tttt_1.addSource('obsExclusionM1','orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ObsM')
T1tttt_1.addSource('obsExclusionP1','orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ObsP')
T1tttt_1.addSource('expExclusion', 'orig/ModelA1.root', 'root', objectName = 'graph_smoothed_Exp')
T1tttt_1.addSource('expExclusionM1', 'orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ExpM')
T1tttt_1.addSource('expExclusionP1', 'orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ExpP')
#----limit source----
T1tttt_1.addSource('efficiencyMap', 'orig/SR23_HighPt.txt', 'txt')
T1tttt_1.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/ModelA1.root'
T1ttttoff.addMassPlane(T1tttt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR24_HighPt")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR24_HighPt", observedN = 7, expectedBG = 2.75 , bgError = 1.18, upperLimit = '5.481E-01*fb', expectedUpperLimit = '2.506E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.constraint = "[[['b','b','W','W']],[['b','b','W','W']]]"
T1ttttoff.conditionDescription ="None"
T1ttttoff.condition = "None"
T1ttttoff.massConstraint = [['dm <= 338.']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure ='Figure 5'
T1tttt_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/Figure5_A1.pdf"
T1tttt_1.addSource('obsExclusion','orig/ModelA1.root', 'root', objectName = 'graph_smoothed_Obs')
T1tttt_1.addSource('obsExclusionM1','orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ObsM')
T1tttt_1.addSource('obsExclusionP1','orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ObsP')
T1tttt_1.addSource('expExclusion', 'orig/ModelA1.root', 'root', objectName = 'graph_smoothed_Exp')
T1tttt_1.addSource('expExclusionM1', 'orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ExpM')
T1tttt_1.addSource('expExclusionP1', 'orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ExpP')
#----limit source----
T1tttt_1.addSource('efficiencyMap', 'orig/SR24_HighPt.txt', 'txt')
T1tttt_1.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/ModelA1.root'
T1ttttoff.addMassPlane(T1tttt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR28_HighPt")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR28_HighPt", observedN = 2, expectedBG = 2.15 , bgError = 0.98, upperLimit = '2.575E-01*fb', expectedUpperLimit = '2.585E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.constraint = "[[['b','b','W','W']],[['b','b','W','W']]]"
T1ttttoff.conditionDescription ="None"
T1ttttoff.condition = "None"
T1ttttoff.massConstraint = [['dm <= 338.']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure ='Figure 5'
T1tttt_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/T1tttt_mostSensitive.pdf"
T1tttt_1.addSource('obsExclusion','orig/ModelA1.root', 'root', objectName = 'graph_smoothed_Obs')
T1tttt_1.addSource('obsExclusionM1','orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ObsM')
T1tttt_1.addSource('obsExclusionP1','orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ObsP')
T1tttt_1.addSource('expExclusion', 'orig/ModelA1.root', 'root', objectName = 'graph_smoothed_Exp')
T1tttt_1.addSource('expExclusionM1', 'orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ExpM')
T1tttt_1.addSource('expExclusionP1', 'orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ExpP')
#----limit source----
T1tttt_1.addSource('efficiencyMap', 'orig/SR28_HighPt.txt', 'txt')
T1tttt_1.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/ModelA1.root'
T1ttttoff.addMassPlane(T1tttt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR21_HighPt")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR21_HighPt", observedN = 12, expectedBG = 7.06 , bgError = 2.5, upperLimit = '7.015E-01*fb', expectedUpperLimit = '4.335E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.constraint = "[[['b','b','W','W']],[['b','b','W','W']]]"
T1ttttoff.conditionDescription ="None"
T1ttttoff.condition = "None"
T1ttttoff.massConstraint = [['dm <= 338.']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure ='Figure 5'
T1tttt_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/Figure5_A1.pdf"
T1tttt_1.addSource('obsExclusion','orig/ModelA1.root', 'root', objectName = 'graph_smoothed_Obs')
T1tttt_1.addSource('obsExclusionM1','orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ObsM')
T1tttt_1.addSource('obsExclusionP1','orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ObsP')
T1tttt_1.addSource('expExclusion', 'orig/ModelA1.root', 'root', objectName = 'graph_smoothed_Exp')
T1tttt_1.addSource('expExclusionM1', 'orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ExpM')
T1tttt_1.addSource('expExclusionP1', 'orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ExpP')
#----limit source----
T1tttt_1.addSource('efficiencyMap', 'orig/SR21_HighPt.txt', 'txt')
T1tttt_1.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/ModelA1.root'
T1ttttoff.addMassPlane(T1tttt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR22_HighPt")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR22_HighPt", observedN = 1, expectedBG = 0.96 , bgError = 0.55, upperLimit = '2.154E-01*fb', expectedUpperLimit = '1.549E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.massConstraint = None
T1tttt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.constraint = "[[['b','b','W','W']],[['b','b','W','W']]]"
T1ttttoff.conditionDescription ="None"
T1ttttoff.condition = "None"
T1ttttoff.massConstraint = [['dm <= 338.']]*2
T1ttttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure ='Figure 5'
T1tttt_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/Figure5_A1.pdf"
T1tttt_1.addSource('obsExclusion','orig/ModelA1.root', 'root', objectName = 'graph_smoothed_Obs')
T1tttt_1.addSource('obsExclusionM1','orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ObsM')
T1tttt_1.addSource('obsExclusionP1','orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ObsP')
T1tttt_1.addSource('expExclusion', 'orig/ModelA1.root', 'root', objectName = 'graph_smoothed_Exp')
T1tttt_1.addSource('expExclusionM1', 'orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ExpM')
T1tttt_1.addSource('expExclusionP1', 'orig/ModelA1.root', 'root', objectName = 'graph_smoothed_ExpP')
#----limit source----
T1tttt_1.addSource('efficiencyMap', 'orig/SR22_HighPt.txt', 'txt')
T1tttt_1.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13013/ModelA1.root'
T1ttttoff.addMassPlane(T1tttt_1)

databaseCreator.create()
