#!/usr/bin/env python3

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
info = MetaInfoInput('CMS-SUS-14-021')
info.url = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS14021'
info.sqrts = 8
info.lumi = 19.7
info.prettyName = 'soft leptons, low jet multiplicity, high ETmiss'
info.private = False
info.arxiv = 'http://arxiv.org/abs/1512.08002'
info.comment ='stop 4 body decay'




#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SRSL1c")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SRSL1c", observedN = 12, expectedBG = 12.3 , bgError = 4.0, upperLimit = '5.787E-01*fb', expectedUpperLimit = '5.777E-01*fb')
#+++++++ next txName block ++++++++++++++
T2bbWW =  dataset.addTxName('T2bbWW')
T2bbWW.constraint = "[[['b','W']],[['b','W']]]"
T2bbWW.conditionDescription = None
T2bbWW.condition =None
T2bbWW.massConstraint = None
T2bbWW.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T2bbWWoff = dataset.addTxName('T2bbWWoff')
T2bbWWoff.constraint = "6.94*[[['b','mu','nu']],[['b','jet','jet']]]"
T2bbWWoff.conditionDescription =None
T2bbWWoff.condition =None
T2bbWWoff.massConstraint = [['dm <= 76.']]*2
T2bbWWoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T2bbWW_1 = T2bbWW.addMassPlane([[x,y]]*2)
#---- new efficiency map -----
#----figure----
T2bbWW_1.figure  = "efficienciesSRSL1c"
T2bbWW_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS14021/efficienciesSRSL1c.png"
T2bbWW_1.addSource('obsExclusion','orig/singleLeptonLimitHistograms.root', 'root', objectName = 'gObserved')
T2bbWW_1.addSource('obsExclusionM1','orig/singleLeptonLimitHistograms.root', 'root', objectName = 'gObservedDown')
T2bbWW_1.addSource('obsExclusionP1','orig/singleLeptonLimitHistograms.root', 'root', objectName = 'gObservedUp')
T2bbWW_1.addSource('expExclusion','orig/singleLeptonLimitHistograms.root', 'root', objectName = 'gExpected')
T2bbWW_1.addSource('expExclusionM1','orig/singleLeptonLimitHistograms.root', 'root', objectName = 'gExpectedDown')
T2bbWW_1.addSource('expExclusionP1','orig/singleLeptonLimitHistograms.root', 'root', objectName = 'gExpectedUp')
T2bbWW_1.addSource('efficiencyMap','orig/efficienciesSRSL.root', 'root', objectName = 'effSRSL1c',scale = 0.01)
T2bbWW_1.dataUrl  = "https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS14021"
T2bbWWoff.addMassPlane(T2bbWW_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SRSL1b")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SRSL1b", observedN = 59, expectedBG = 81.3 , bgError = 19.1, upperLimit = '1.473E+00*fb', expectedUpperLimit = '2.148E+00*fb')
#+++++++ next txName block ++++++++++++++
T2bbWW =  dataset.addTxName('T2bbWW')
T2bbWW.constraint = "[[['b','W']],[['b','W']]]"
T2bbWW.conditionDescription = None
T2bbWW.condition =None
T2bbWW.massConstraint = None
T2bbWW.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T2bbWWoff = dataset.addTxName('T2bbWWoff')
T2bbWWoff.constraint = "6.94*[[['b','mu','nu']],[['b','jet','jet']]]"
T2bbWWoff.conditionDescription =None
T2bbWWoff.condition =None
T2bbWWoff.massConstraint = [['dm <= 76.']]*2
T2bbWWoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T2bbWW_1 = T2bbWW.addMassPlane([[x,y]]*2)
#---- new efficiency map -----
#----figure----
T2bbWW_1.figure  = "efficienciesSRSL1b"
T2bbWW_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS14021/efficienciesSRSL1b.png"
T2bbWW_1.addSource('obsExclusion','orig/singleLeptonLimitHistograms.root', 'root', objectName = 'gObserved')
T2bbWW_1.addSource('obsExclusionM1','orig/singleLeptonLimitHistograms.root', 'root', objectName = 'gObservedDown')
T2bbWW_1.addSource('obsExclusionP1','orig/singleLeptonLimitHistograms.root', 'root', objectName = 'gObservedUp')
T2bbWW_1.addSource('expExclusion','orig/singleLeptonLimitHistograms.root', 'root', objectName = 'gExpected')
T2bbWW_1.addSource('expExclusionM1','orig/singleLeptonLimitHistograms.root', 'root', objectName = 'gExpectedDown')
T2bbWW_1.addSource('expExclusionP1','orig/singleLeptonLimitHistograms.root', 'root', objectName = 'gExpectedUp')
T2bbWW_1.addSource('efficiencyMap','orig/efficienciesSRSL.root', 'root', objectName = 'effSRSL1b',scale = 0.01)
T2bbWW_1.dataUrl  = "https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS14021"
T2bbWWoff.addMassPlane(T2bbWW_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SRSL1a")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SRSL1a", observedN = 121, expectedBG = 134.5 , bgError = 19.8, upperLimit = '1.880E+00*fb', expectedUpperLimit = '2.318E+00*fb')
#+++++++ next txName block ++++++++++++++
T2bbWW =  dataset.addTxName('T2bbWW')
T2bbWW.constraint = "[[['b','W']],[['b','W']]]"
T2bbWW.conditionDescription = None
T2bbWW.condition =None
T2bbWW.massConstraint = None
T2bbWW.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T2bbWWoff = dataset.addTxName('T2bbWWoff')
T2bbWWoff.constraint = "6.94*[[['b','mu','nu']],[['b','jet','jet']]]"
T2bbWWoff.conditionDescription =None
T2bbWWoff.condition =None
T2bbWWoff.massConstraint = [['dm <= 76.']]*2
T2bbWWoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T2bbWW_1 = T2bbWW.addMassPlane([[x,y]]*2)
#---- new efficiency map -----
#----figure----
T2bbWW_1.figure  = "efficienciesSRSL1a"
T2bbWW_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS14021/efficienciesSRSL1a.png"
T2bbWW_1.addSource('obsExclusion','orig/singleLeptonLimitHistograms.root', 'root', objectName = 'gObserved')
T2bbWW_1.addSource('obsExclusionM1','orig/singleLeptonLimitHistograms.root', 'root', objectName = 'gObservedDown')
T2bbWW_1.addSource('obsExclusionP1','orig/singleLeptonLimitHistograms.root', 'root', objectName = 'gObservedUp')
T2bbWW_1.addSource('expExclusion','orig/singleLeptonLimitHistograms.root', 'root', objectName = 'gExpected')
T2bbWW_1.addSource('expExclusionM1','orig/singleLeptonLimitHistograms.root', 'root', objectName = 'gExpectedDown')
T2bbWW_1.addSource('expExclusionP1','orig/singleLeptonLimitHistograms.root', 'root', objectName = 'gExpectedUp')
T2bbWW_1.addSource('efficiencyMap','orig/efficienciesSRSL.root', 'root', objectName = 'effSRSL1a',scale = 0.01)
T2bbWW_1.dataUrl  = "https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS14021"
T2bbWWoff.addMassPlane(T2bbWW_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SRSL2")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SRSL2", observedN = 41, expectedBG = 32.1 , bgError = 7.7, upperLimit = '1.400E+00*fb', expectedUpperLimit = '1.010E+00*fb')
#+++++++ next txName block ++++++++++++++
T2bbWW =  dataset.addTxName('T2bbWW')
T2bbWW.constraint = "[[['b','W']],[['b','W']]]"
T2bbWW.conditionDescription = None
T2bbWW.condition =None
T2bbWW.massConstraint = None
T2bbWW.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T2bbWWoff = dataset.addTxName('T2bbWWoff')
T2bbWWoff.constraint = "6.94*[[['b','mu','nu']],[['b','jet','jet']]]"
T2bbWWoff.conditionDescription =None
T2bbWWoff.condition =None
T2bbWWoff.massConstraint = [['dm <= 76.']]*2
T2bbWWoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T2bbWW_1 = T2bbWW.addMassPlane([[x,y]]*2)
#---- new efficiency map -----
#----figure----
T2bbWW_1.figure  = "efficienciesSRSL2"
T2bbWW_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS14021/efficienciesSRSL2.png"
T2bbWW_1.addSource('obsExclusion','orig/singleLeptonLimitHistograms.root', 'root', objectName = 'gObserved')
T2bbWW_1.addSource('obsExclusionM1','orig/singleLeptonLimitHistograms.root', 'root', objectName = 'gObservedDown')
T2bbWW_1.addSource('obsExclusionP1','orig/singleLeptonLimitHistograms.root', 'root', objectName = 'gObservedUp')
T2bbWW_1.addSource('expExclusion','orig/singleLeptonLimitHistograms.root', 'root', objectName = 'gExpected')
T2bbWW_1.addSource('expExclusionM1','orig/singleLeptonLimitHistograms.root', 'root', objectName = 'gExpectedDown')
T2bbWW_1.addSource('expExclusionP1','orig/singleLeptonLimitHistograms.root', 'root', objectName = 'gExpectedUp')
T2bbWW_1.addSource('efficiencyMap','orig/efficienciesSRSL.root', 'root', objectName = 'effSRSL2',scale = 0.01)
T2bbWW_1.dataUrl  = "https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS14021"
T2bbWWoff.addMassPlane(T2bbWW_1)

databaseCreator.create()
