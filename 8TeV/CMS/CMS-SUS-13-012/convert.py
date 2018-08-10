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
info = MetaInfoInput('CMS-SUS-13-012')
info.sqrts = '8.0'
info.private = False
info.lumi = '19.5'
info.publication = 'JHEP06(2014)055'
info.url = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13012'
info.arxiv = 'http://arxiv.org/abs/1402.4770'
info.contact = 'cms-pag-conveners-sus@NOSPAMcernSPAMNOT.ch'
info.prettyName = 'jet multiplicity + HTmiss'
info.implementedBy = 'Federico A.'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T2 = dataset.addTxName('T2')
T2.constraint ="[[['jet']],[['jet']]]"
T2.conditionDescription ="None"
T2.condition ="None"
T2.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T2 = T2.addMassPlane(2*[[x, y]])
T2.figure = "Fig_7a"
T2.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13012/Fig_7a.pdf"
T2.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13012/SUS13012_XsecLimits_T2qq.root"
T2.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'expectedUpperLimits', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/SUS13012_XsecLimits_T2qq.root', 'orig/SUS13012_XsecLimits_T2qq.root', 'orig/SUS13012_XsecLimits_T2qq.root', 'orig/SUS13012_XsecLimits_T2qq.root', 'orig/SUS13012_XsecLimits_T2qq.root', 'orig/SUS13012_XsecLimits_T2qq.root', 'orig/SUS13012_XsecLimits_T2qq.root', 'orig/SUS13012_XsecLimits_T2qq.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['combined_expExclOneTimesProspino', 'combined_expExclMinusOneSigmaProspino', 'combined_expExclPlusOneSigmaProspino', 'combined_expLimit', 'combined_obsExclOneTimesProspino', 'combined_obsExclMinusSysErrProspino', 'combined_obsExclPlusSysErrProspino', 'combined_obsLimit'])

#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.source = "CMS"
T1tttt.massConstraint = None
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.constraint = "[[['b','W','b','W']],[['b','W','b','W']]]"
T1ttttoff.conditionDescription ="None"
T1ttttoff.condition = "None"
T1ttttoff.massConstraint = [['dm <= 338.0'], ['dm <= 338.0']]
T1ttttoff.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane(2*[[x, y]])
T1tttt.figure = "Fig_7c"
T1tttt.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13012/Fig_7c.pdf"
T1tttt.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13012/SUS13012_XsecLimits_T1tttt.root"
T1tttt.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'expectedUpperLimits', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/SUS13012_XsecLimits_T1tttt.root', 'orig/SUS13012_XsecLimits_T1tttt.root', 'orig/SUS13012_XsecLimits_T1tttt.root', 'orig/SUS13012_XsecLimits_T1tttt.root', 'orig/SUS13012_XsecLimits_T1tttt.root', 'orig/SUS13012_XsecLimits_T1tttt.root', 'orig/SUS13012_XsecLimits_T1tttt.root', 'orig/SUS13012_XsecLimits_T1tttt.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['combined_expExclOneTimesProspino', 'combined_expExclMinusOneSigmaProspino', 'combined_expExclPlusOneSigmaProspino', 'combined_expLimit', 'combined_obsExclOneTimesProspino', 'combined_obsExclMinusSysErrProspino', 'combined_obsExclPlusSysErrProspino', 'combined_obsLimit'])
T1ttttoff.addMassPlane(T1tttt)

#+++++++ next txName block ++++++++++++++
T1 = dataset.addTxName('T1')
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription ="None"
T1.condition ="None"
T1.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T1 = T1.addMassPlane(2*[[x, y]])
T1.figure = "Fig_7b"
T1.figureUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13012/Fig_7b.pdf"
T1.dataUrl = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13012/SUS13012_XsecLimits_T1qqqq.root"
T1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'expectedUpperLimits', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/SUS13012_XsecLimits_T1qqqq.root', 'orig/SUS13012_XsecLimits_T1qqqq.root', 'orig/SUS13012_XsecLimits_T1qqqq.root', 'orig/SUS13012_XsecLimits_T1qqqq.root', 'orig/SUS13012_XsecLimits_T1qqqq.root', 'orig/SUS13012_XsecLimits_T1qqqq.root', 'orig/SUS13012_XsecLimits_T1qqqq.root', 'orig/SUS13012_XsecLimits_T1qqqq.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['combined_expExclOneTimesProspino', 'combined_expExclMinusOneSigmaProspino', 'combined_expExclPlusOneSigmaProspino', 'combined_expLimit', 'combined_obsExclOneTimesProspino', 'combined_obsExclMinusSysErrProspino', 'combined_obsExclPlusSysErrProspino', 'combined_obsLimit'])



databaseCreator.create()
