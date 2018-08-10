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
info = MetaInfoInput('CMS-PAS-SUS-13-018')
info.url = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13018'
info.sqrts = 8
info.lumi = 19.4
info.prettyName = '1-2 b-jets + Etmiss, M_CT'
info.private = False
info.arxiv = ''
info.contact = ''
info.publication = ''
info.comment = 'Only https://cds.cern.ch/record/1693164 document, will be superseeded by CMS-SUS-14-001'
info.supersedes =''
info.supersededBy = ''


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T2bb = dataset.addTxName('T2bb')
T2bb.checked =''
T2bb.constraint = "[[['b']],[['b']]]"
T2bb.conditionDescription = None
T2bb.condition = None
T2bb.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T2bb_1 = T2bb.addMassPlane(2*[[x, y]])
T2bb_1.figure = 'Fig. 6'
T2bb_1.figureUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13018/T2bbXSEC.png'
T2bb_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13018'
T2bb_1.histoDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13018/T2bbXSEC.C'
T2bb_1.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13018/T2bbXSEC.C'
T2bb_1.exclusionDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13018/T2bbXSEC.C'
T2bb_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/combined_expExclOneTimesProspino.root', 'orig/combined_expExclMinusOneSigmaProspino.root', 'orig/combined_expExclPlusOneSigmaProspino.root', 'orig/combined_obsExclOneTimesProspino.root', 'orig/combined_obsExclMinusSysErrProspino.root', 'orig/combined_obsExclPlusSysErrProspino.root', 'orig/combined_expLimit.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['combined_expExclOneTimesProspino', 'combined_expExclMinusOneSigmaProspino', 'combined_expExclPlusOneSigmaProspino', 'combined_obsExclOneTimesProspino', 'combined_obsExclMinusSysErrProspino', 'combined_obsExclPlusSysErrProspino', 'combined_expLimit'])



databaseCreator.create()
