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
info = MetaInfoInput('CMS-SUS-13-002')
info.url = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13002'
info.sqrts = 8
info.lumi =  19.5
info.prettyName = '>= 3 leptons (+jets) + Etmiss'
info.private = False
info.arxiv = 'http://arxiv.org/abs/1404.5801v2'
info.contact =''
info.publication = 'http://journals.aps.org/prd/abstract/10.1103/PhysRevD.90.032006'
info.comment =''
info.supersedes = 'CMS-PAS-SUS-12-026'
info.supersededBy =''


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked =''
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = None
T1tttt.condition =None
T1tttt.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane(2*[[x, y]])
T1tttt_1.figure = 'Figure 11'
T1tttt_1.figureUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13002/curve_T1tttt_overlay_observed.pdf'
T1tttt_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13002'
T1tttt_1.histoDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13002/curve_T1tttt_overlay_observed.pdf'
T1tttt_1.dataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13002/Limit_T1tttt.root'
T1tttt_1.exclusionDataUrl = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13002'
T1tttt_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'expectedUpperLimits', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/contours_T1tttt.root', 'orig/contours_T1tttt.root', 'orig/contours_T1tttt.root', 'orig/Limit_T1tttt.root', 'orig/contours_T1tttt.root', 'orig/contours_T1tttt.root', 'orig/contours_T1tttt.root', 'orig/Limit_T1tttt.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['Expected', 'Expected1m', 'Expected1p', 'hrExp', 'Observed', 'Observed1m', 'Observed1p', 'hrObs'],units= [None, None, None, None, None, None, None, 'fb'])



databaseCreator.create()
