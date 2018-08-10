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
info = MetaInfoInput('CMS-PAS-SUS-14-011')
info.url = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS14011'
info.sqrts = 8
info.lumi = 19.3
info.prettyName = 'razor with b-jets'
info.private = False
info.arxiv = ""
info.contact = ""
info.publication = ""
info.comment = 'Only PAS document https://cms-physics.web.cern.ch/cms-physics/public/SUS-14-011-pas.pdf'
info.supersededBy = 'CMS-SUS-13-004'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked =''
T1bbbb.constraint = "[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription = None
T1bbbb.condition = None
T1bbbb.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T1bbbb_1 = T1bbbb.addMassPlane(2*[[x, y]])
T1bbbb_1.figure = 'Fig. 4'
T1bbbb_1.figureUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS14011/T1bbbbHybridNew0LXSEC.pdf'
T1bbbb_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS14011'
T1bbbb_1.histoDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS14011/GluinoLimitsRazor2014.root'
T1bbbb_1.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS14011/GluinoLimitsRazor2014.root'
T1bbbb_1.exclusionDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS14011/GluinoLimitsRazor2014.root'
T1bbbb_1.setSources(dataLabels= ['expExclusion', 'expExclusionP1', 'expectedUpperLimits', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/GluinoLimitsRazor2014.root', 'orig/GluinoLimitsRazor2014.root', 'orig/GluinoLimitsRazor2014.root', 'orig/GluinoLimitsRazor2014.root', 'orig/GluinoLimitsRazor2014.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root'],objectNames= ['Exp_contour_T1bbbb_0L', 'ExpMinus_contour_T1bbbb_0L', 'Exp_xsec_T1bbbb_0L', 'Obs_contour_T1bbbb_0L', 'Obs_xsec_T1bbbb_0L'])

#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.checked =''
T2tt.constraint = "[[['t']],[['t']]]"
T2tt.conditionDescription = None
T2tt.condition = None
T2tt.source = "CMS"
T2tt.massConstraint = None
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.checked =''
T2ttoff.constraint =  "[[['b','W']],[['b','W']]]"
T2ttoff.conditionDescription = None
T2ttoff.condition = None
T2ttoff.massConstraint = [['dm <= 169.0'], ['dm <= 169.0']]
T2ttoff.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane(2*[[x, y]])
T2tt_1.figure = 'Fig. 5'
T2tt_1.figureUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS14011/T2ttHybridNew0Lp1Lp2LXSEC.pdf'
T2tt_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS14011'
T2tt_1.histoDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS14011/StopLimitsRazor2014.root'
T2tt_1.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS14011/StopLimitsRazor2014.root'
T2tt_1.exclusionDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS14011/StopLimitsRazor2014.root'
T2tt_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'expectedUpperLimits', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/StopLimitsRazor2014.root', 'orig/StopLimitsRazor2014.root', 'orig/StopLimitsRazor2014.root', 'orig/StopLimitsRazor2014.root', 'orig/StopLimitsRazor2014.root', 'orig/StopLimitsRazor2014.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['Exp_contour_T2tt_0Lp1Lp2L', 'ExpMinus_contour_T2tt_0Lp1Lp2L', 'ExpPlus_contour_T2tt_0Lp1Lp2L', 'Exp_xsec_T2tt_0Lp1Lp2L', 'Obs_contour_T2tt_0Lp1Lp2L', 'Obs_xsec_T2tt_0Lp1Lp2L'])
T2ttoff.addMassPlane(T2tt_1)

#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked =''
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = None
T1tttt.condition = None
T1tttt.source = "CMS"
T1tttt.massConstraint = None
T1ttttoff = dataset.addTxName('T1ttttoff')
T1ttttoff.checked =''
T1ttttoff.constraint = "[[['b','W','b','W']],[['b','W','b','W']]]"
T1ttttoff.conditionDescription =None
T1ttttoff.condition =None
T1ttttoff.massConstraint = [['dm <= 338.0'], ['dm <= 338.0']]
T1ttttoff.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane(2*[[x, y]])
T1tttt_1.figure = 'Fig. 4'
T1tttt_1.figureUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS14011/T1ttttHybridNew0Lp1Lp2LXSEC.pdf'
T1tttt_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS14011'
T1tttt_1.histoDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS14011/GluinoLimitsRazor2014.root'
T1tttt_1.dataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS14011/GluinoLimitsRazor2014.root'
T1tttt_1.exclusionDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS14011/GluinoLimitsRazor2014.root'
T1tttt_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'expectedUpperLimits', 'obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/GluinoLimitsRazor2014.root', 'orig/GluinoLimitsRazor2014.root', 'orig/GluinoLimitsRazor2014.root', 'orig/GluinoLimitsRazor2014.root', 'orig/GluinoLimitsRazor2014.root', 'orig/GluinoLimitsRazor2014.root'],
                 dataFormats= ['root', 'root', 'root', 'root', 'root', 'root'],objectNames= ['Exp_contour_T1tttt_0Lp1Lp2L', 'ExpMinus_contour_T1tttt_0Lp1Lp2L', 'ExpPlus_contour_T1tttt_0Lp1Lp2L', 'Exp_xsec_T1tttt_0Lp1Lp2L', 'Obs_contour_T1tttt_0Lp1Lp2L', 'Obs_xsec_T1tttt_0Lp1Lp2L'])
T1ttttoff.addMassPlane(T1tttt_1)



databaseCreator.create()
