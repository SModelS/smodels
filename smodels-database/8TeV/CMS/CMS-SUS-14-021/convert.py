#!/usr/bin/env python3

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
info = MetaInfoInput('CMS-SUS-14-021')
info.url = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS14021'
info.sqrts = 8
info.lumi = 19.7
info.prettyName = 'soft leptons, low jet multiplicity, high ETmiss'
info.private = False
info.arxiv = 'http://arxiv.org/abs/1512.08002'
info.contact = ''
info.publication =''
info.comment ='ISR jet. Only single lepton analysis data available and implemented.'
info.supersedes =''
info.supersededBy =''


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T2bbWW = dataset.addTxName('T2bbWW')
T2bbWW.checked = ''
T2bbWW.constraint = "[[['b','W']],[['b','W']]]"
T2bbWW.conditionDescription = None
T2bbWW.condition =None
T2bbWW.source = "CMS"
T2bbWW.massConstraint = None
T2bbWWoff = dataset.addTxName('T2bbWWoff')
T2bbWWoff.checked = ''
T2bbWWoff.constraint = "3.47*[[['b','l','nu']],[['b','jet','jet']]]"
T2bbWWoff.conditionDescription ="None"
T2bbWWoff.condition ="None"
T2bbWWoff.massConstraint = [['dm <= 76.0'], ['dm <= 76.0']]
T2bbWWoff.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T2bbWW_1 = T2bbWW.addMassPlane(2*[[x, y]])
T2bbWW_1.figure = 'Fig.3 Left'
T2bbWW_1.figureUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS14021/singleLeptonLimits.pdf'
T2bbWW_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS14021'
T2bbWW_1.histoDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS14021/singleLeptonLimitHistograms.root'
T2bbWW_1.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS14021/singleLeptonLimitHistograms.root'
T2bbWW_1.exclusionDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS14021/singleLeptonLimitHistograms.root'
T2bbWW_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'expectedUpperLimits', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/singleLeptonLimitHistograms.root', 'orig/singleLeptonLimitHistograms.root', 'orig/singleLeptonLimitHistograms.root', 'orig/singleLeptonLimitHistograms.root', 'orig/singleLeptonLimitHistograms.root', 'orig/singleLeptonLimitHistograms.root', 'orig/singleLeptonLimitHistograms.root', 'orig/singleLeptonLimitHistograms.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['gExpected', 'gExpectedDown', 'gExpectedUp', 'hExpected', 'gObserved', 'gObservedDown', 'gObservedUp', 'hObserved'])
T2bbWWoff.addMassPlane(T2bbWW_1)



databaseCreator.create()
