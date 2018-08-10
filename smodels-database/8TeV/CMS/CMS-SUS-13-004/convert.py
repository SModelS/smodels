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
info = MetaInfoInput('CMS-SUS-13-004')
info.url = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13004'
info.sqrts = 8
info.lumi = 19.3
info.prettyName = '>= 1 b-jet + Etmiss, Razor'
info.private = False
info.arxiv = 'http://arxiv.org/abs/1502.00300v2'
info.publication = 'http://dx.doi.org/10.1103/PhysRevD.91.052018'
info.comment ='Other topologies have a compressed mass spectrum (decays with intermediate charginos, <= 5 GeV mass difference)'
info.supersedes = 'CMS-PAS-SUS-14-011'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

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
T2tt_1.figure = 'Fig. 15c'
T2tt_1.figureUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13004/T2ttCOMBINED.pdf'
T2tt_1.dataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13004/razor_8TeV_sms_results.root'
T2tt_1.histoDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13004/razor_8TeV_sms_results.root'
T2tt_1.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13004/razor_8TeV_sms_results.root'
T2tt_1.exclusionDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13004/razor_8TeV_sms_results.root'
T2tt_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/T2tt_Comb.root', 'orig/T2tt_Comb.root', 'orig/T2tt_Comb.root', 'orig/T2tt_Comb.root', 'orig/T2tt_Comb.root', 'orig/T2tt_Comb.root', 'orig/T2tt_Comb.root'],
                 dataFormats= ['canvas', 'canvas', 'canvas', 'canvas', 'canvas', 'canvas', 'canvas'],objectNames= ['cCONT_', 'cCONT_', 'cCONT_', 'cCONT_', 'cCONT_', 'cCONT_', 'cCONT_'],indices= [5, 7, 6, 8, 10, 9, 2])
T2ttoff.addMassPlane(T2tt_1)

#+++++++ next txName block ++++++++++++++
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked =''
T1bbbb.constraint = "[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription = None
T1bbbb.condition = None
T1bbbb.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
T1bbbb_1 = T1bbbb.addMassPlane(2*[[x, y]])
T1bbbb_1.figure = 'Fig. 13a'
T1bbbb_1.figureUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13004/T1bbbbHybridNew0LXSEC.png'
T1bbbb_1.dataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13004/razor_8TeV_sms_results.root'
T1bbbb_1.histoDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13004/razor_8TeV_sms_results.root'
T1bbbb_1.dataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13004/razor_8TeV_sms_results.root'
T1bbbb_1.exclusionDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13004/razor_8TeV_sms_results.root'
T1bbbb_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/T1bbbb.root', 'orig/T1bbbb.root', 'orig/T1bbbb.root', 'orig/T1bbbb.root', 'orig/T1bbbb.root', 'orig/T1bbbb.root', 'orig/T1bbbb.root'],
                 dataFormats= ['canvas', 'canvas', 'canvas', 'canvas', 'canvas', 'canvas', 'canvas'],objectNames= ['cCONT_', 'cCONT_', 'cCONT_', 'cCONT_', 'cCONT_', 'cCONT_', 'cCONT_'],indices= [5, 7, 6, 8, 10, 9, 2])

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
T1tttt_1.figure = 'Fig. 13e'
T1tttt_1.figureUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13004/T1ttttHybridNew0Lp1Lp2LXSEC.pdf'
T1tttt_1.dataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13004/razor_8TeV_sms_results.root'
T1tttt_1.histoDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13004/razor_8TeV_sms_results.root'
T1tttt_1.dataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13004/razor_8TeV_sms_results.root'
T1tttt_1.exclusionDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13004/razor_8TeV_sms_results.root'
T1tttt_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/T1tttt.root', 'orig/T1tttt.root', 'orig/T1tttt.root', 'orig/T1tttt.root', 'orig/T1tttt.root', 'orig/T1tttt.root', 'orig/T1tttt.root'],
                 dataFormats= ['canvas', 'canvas', 'canvas', 'canvas', 'canvas', 'canvas', 'canvas'],objectNames= ['cCONT_', 'cCONT_', 'cCONT_', 'cCONT_', 'cCONT_', 'cCONT_', 'cCONT_'],indices= [5, 7, 6, 8, 10, 9, 2])
T1ttttoff.addMassPlane(T1tttt_1)



databaseCreator.create()
