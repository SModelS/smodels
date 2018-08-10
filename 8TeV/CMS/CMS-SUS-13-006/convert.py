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
argparser.add_argument ('-t', '--ntoys',
    help = 'number of toys to throw',\
    type = int, default=200000  )
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

DataSetInput.ntoys = args.ntoys


#+++++++ global info block ++++++++++++++
info = MetaInfoInput('CMS-SUS-13-006')
info.url ='https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13006'
info.sqrts = 8
info.lumi = 19.5
info.prettyName = 'EW productions with decays to leptons, W, Z, and Higgs'
info.private = False
info.arxiv = 'http://arxiv.org/abs/1405.7570'
info.contact = ''
info.publication = 'http://link.springer.com/article/10.1140%2Fepjc%2Fs10052-014-3036-7'
info.comment ='TSlepSlep only LeftHanded'
info.supersedes ='CMS-PAS-SUS-12-022'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
TChiWZ = dataset.addTxName('TChiWZ')
TChiWZ.checked =''
TChiWZ.constraint = "[[['W']],[['Z']]]"
TChiWZ.conditionDescription =None
TChiWZ.condition =None
TChiWZ.source = "CMS"
TChiWZ.massConstraint = None
TChiWZoff = dataset.addTxName('TChiWZoff')
TChiWZoff.checked =''
TChiWZoff.constraint = "71.*([[['mu+','mu-']],[['l','nu']]] + [[['e+','e-']],[['l','nu']]])"
TChiWZoff.conditionDescription =None
TChiWZoff.condition = "cGtr([[['mu+','mu-']],[['l','nu']]],[[['e+','e-']],[['l','nu']]])"
TChiWZoff.massConstraint = [['dm <= 86.0'], ['dm <= 86.0']]
TChiWZoff.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
TChiWZ_1 = TChiWZ.addMassPlane(2*[[x, y]])
TChiWZ_1.figure = 'Fig. 16(left)'
TChiWZ_1.figureUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/Fig16_exclusion_TChiWZ.png'
TChiWZ_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13006'
TChiWZ_1.histoDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/exclusion_TChiWZ.root'
TChiWZ_1.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/exclusion_TChiWZ.root'
TChiWZ_1.exclusionDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/exclusion_TChiWZ.root'
TChiWZ_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/exclusion_TChiWZ.root', 'orig/exclusion_TChiWZ.root', 'orig/exclusion_TChiWZ.root', 'orig/exclusion_TChiWZ.root', 'orig/exclusion_TChiWZ.root', 'orig/exclusion_TChiWZ.root', 'orig/exclusion_TChiWZ.root'],
                 dataFormats= ['canvas', 'canvas', 'canvas', 'canvas', 'canvas', 'canvas', 'canvas'],objectNames= ['interpret', 'interpret', 'interpret', 'interpret', 'interpret', 'interpret', 'interpret'],indices= [9, 5, 6, 12, 8, 7, 2],units= [None, None, None, None, None, None, 'fb'])
TChiWZoff.addMassPlane(TChiWZ_1)

#+++++++ next txName block ++++++++++++++
TChiChipmSlepStau = dataset.addTxName('TChiChipmSlepStau')
TChiChipmSlepStau.checked =''
TChiChipmSlepStau.constraint = "[[['L'],['L']],[['nu'],['ta']]]"
TChiChipmSlepStau.conditionDescription = None
TChiChipmSlepStau.condition =  "cGtr([[['L'],['L']],[['nu'],['ta']]],3.*[[['ta'],['ta']],[['nu'],['ta']]]);cGtr([[['L'],['L']],[['nu'],['ta']]],3.*[[['e'],['e']],[['nu'],['ta']]])"
TChiChipmSlepStau.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
TChiChipmSlepStau_1 = TChiChipmSlepStau.addMassPlane(2*[[x, 0.05*x+0.95*y, y]])
TChiChipmSlepStau_1.figure ='Fig. 14(top-left)'
TChiChipmSlepStau_1.figureUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/Fig14_exclusion_TChiSlepSnu_2a_0_05.pdf'
TChiChipmSlepStau_1.dataUrl ='https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13006'
TChiChipmSlepStau_1.histoDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/exclusion_TChiSlepSnu_2a_0_05.root'
TChiChipmSlepStau_1.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/exclusion_TChiSlepSnu_2a_0_05.root'
TChiChipmSlepStau_1.exclusionDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/exclusion_TChiSlepSnu_2a_0_05.root'
TChiChipmSlepStau_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/exclusion_TChiSlepSnu_2a_0_05.root', 'orig/exclusion_TChiSlepSnu_2a_0_05.root', 'orig/exclusion_TChiSlepSnu_2a_0_05.root', 'orig/exclusion_TChiSlepSnu_2a_0_05.root', 'orig/exclusion_TChiSlepSnu_2a_0_05.root', 'orig/exclusion_TChiSlepSnu_2a_0_05.root', 'orig/exclusion_TChiSlepSnu_2a_0_05.root'],
                 dataFormats= ['canvas', 'canvas', 'canvas', 'canvas', 'canvas', 'canvas', 'canvas'],objectNames= ['interpret', 'interpret', 'interpret', 'interpret', 'interpret', 'interpret', 'interpret'],indices= [8, 4, 5, 11, 7, 6, 2],units= [None, None, None, None, None, None, 'fb'])
#+++++++ next mass plane block ++++++++++++++
TChiChipmSlepStau_2 = TChiChipmSlepStau.addMassPlane(2*[[x, 0.5*x+0.5*y, y]])
TChiChipmSlepStau_2.figure ='Fig. 14(top-right)'
TChiChipmSlepStau_2.figureUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/Fig14_exclusion_TChiSlepSnu_2a_0_5.pdf'
TChiChipmSlepStau_2.dataUrl ='https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13006'
TChiChipmSlepStau_2.histoDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/exclusion_TChiSlepSnu_2a_0_5.root'
TChiChipmSlepStau_2.dataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/exclusion_TChiSlepSnu_2a_0_5.root'
TChiChipmSlepStau_2.exclusionDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/exclusion_TChiSlepSnu_2a_0_5.root'
TChiChipmSlepStau_2.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/exclusion_TChiSlepSnu_2a_0_5.root', 'orig/exclusion_TChiSlepSnu_2a_0_5.root', 'orig/exclusion_TChiSlepSnu_2a_0_5.root', 'orig/exclusion_TChiSlepSnu_2a_0_5.root', 'orig/exclusion_TChiSlepSnu_2a_0_5.root', 'orig/exclusion_TChiSlepSnu_2a_0_5.root', 'orig/exclusion_TChiSlepSnu_2a_0_5.root'],
                 dataFormats= ['canvas', 'canvas', 'canvas', 'canvas', 'canvas', 'canvas', 'canvas'],objectNames= ['interpret', 'interpret', 'interpret', 'interpret', 'interpret', 'interpret', 'interpret'],indices= [8, 4, 5, 9, 7, 6, 2],units= [None, None, None, None, None, None, 'fb'])
#+++++++ next mass plane block ++++++++++++++
TChiChipmSlepStau_3 = TChiChipmSlepStau.addMassPlane(2*[[x, 0.95*x+0.05*y, y]])
TChiChipmSlepStau_3.figure ='Fig. 14(bottom)'
TChiChipmSlepStau_3.figureUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/Fig14_exclusion_TChiSlepSnu_2a_0_95.pdf'
TChiChipmSlepStau_3.dataUrl ='https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13006'
TChiChipmSlepStau_3.histoDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/exclusion_TChiSlepSnu_2a_0_95.root'
TChiChipmSlepStau_3.dataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/exclusion_TChiSlepSnu_2a_0_95.root'
TChiChipmSlepStau_3.exclusionDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/exclusion_TChiSlepSnu_2a_0_95.root'
TChiChipmSlepStau_3.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/exclusion_TChiSlepSnu_2a_0_95.root', 'orig/exclusion_TChiSlepSnu_2a_0_95.root', 'orig/exclusion_TChiSlepSnu_2a_0_95.root', 'orig/exclusion_TChiSlepSnu_2a_0_95.root', 'orig/exclusion_TChiSlepSnu_2a_0_95.root', 'orig/exclusion_TChiSlepSnu_2a_0_95.root', 'orig/exclusion_TChiSlepSnu_2a_0_95.root'],
                 dataFormats= ['canvas', 'canvas', 'canvas', 'canvas', 'canvas', 'canvas', 'canvas'],objectNames= ['interpret', 'interpret', 'interpret', 'interpret', 'interpret', 'interpret', 'interpret'],indices= [8, 4, 5, 11, 7, 6, 2],units= [None, None, None, None, None, None, 'fb'])

#+++++++ next txName block ++++++++++++++
TSlepSlep = dataset.addTxName('TSlepSlep')
TSlepSlep.checked =''
TSlepSlep.constraint ="[[['e+']],[['e-']]]+[[['mu+']],[['mu-']]]"
TSlepSlep.conditionDescription ="[[['mu+']],[['mu-']]] > [[['e+']],[['e-']]]"
TSlepSlep.condition ="Cgtr([[['mu+']],[['mu-']]],[[['e+']],[['e-']]])"
TSlepSlep.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
TSlepSlep_1 = TSlepSlep.addMassPlane(2*[[x, y]])
TSlepSlep_1.figure ='Fig.18 (upper-right)'
TSlepSlep_1.figureUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/Fig18_exclusion_TSlepSlepL.png'
TSlepSlep_1.dataUrl ='https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13006'
TSlepSlep_1.histoDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/exclusion_TSlepSlepL.root'
TSlepSlep_1.dataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/exclusion_TSlepSlepL.root'
TSlepSlep_1.exclusionDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/exclusion_TSlepSlepL.root'
TSlepSlep_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/exclusion_TSlepSlepL.root', 'orig/exclusion_TSlepSlepL.root', 'orig/exclusion_TSlepSlepL.root', 'orig/exclusion_TSlepSlepL.root', 'orig/exclusion_TSlepSlepL.root', 'orig/exclusion_TSlepSlepL.root', 'orig/exclusion_TSlepSlepL.root'],
                 dataFormats= ['canvas', 'canvas', 'canvas', 'canvas', 'canvas', 'canvas', 'canvas'],objectNames= ['interpret', 'interpret', 'interpret', 'interpret', 'interpret', 'interpret', 'interpret'],indices= [7, 3, 4, 8, 6, 5, 2],units= [None, None, None, None, None, None, 'fb'])

#+++++++ next txName block ++++++++++++++
TChiWH = dataset.addTxName('TChiWH')
TChiWH.checked =''
TChiWH.constraint = "[[['W']],[['higgs']]]"
TChiWH.conditionDescription =None
TChiWH.condition =None
TChiWH.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
TChiWH_1 = TChiWH.addMassPlane(2*[[x, y]])
TChiWH_1.figure = 'Fig. 16(right)'
TChiWH_1.figureUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/Fig16_exclusion_TChiWH.png'
TChiWH_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13006'
TChiWH_1.histoDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/exclusion_TChiWH.root'
TChiWH_1.dataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/exclusion_TChiWH.root'
TChiWH_1.exclusionDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/exclusion_TChiWH.root'
TChiWH_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/exclusion_TChiWH.root', 'orig/exclusion_TChiWH.root', 'orig/exclusion_TChiWH.root', 'orig/exclusion_TChiWH.root', 'orig/exclusion_TChiWH.root', 'orig/exclusion_TChiWH.root', 'orig/exclusion_TChiWH.root'],
                 dataFormats= ['canvas', 'canvas', 'canvas', 'canvas', 'canvas', 'canvas', 'canvas'],objectNames= ['interpret', 'interpret', 'interpret', 'interpret', 'interpret', 'interpret', 'interpret'],indices= [8, 4, 5, 9, 7, 6, 2],units= [None, None, None, None, None, None, 'fb'])

#+++++++ next txName block ++++++++++++++
TChiChipmSlepL = dataset.addTxName('TChiChipmSlepL')
TChiChipmSlepL.checked =''
TChiChipmSlepL.constraint = "2.*([[['L'],['L']],[['L'],['nu']]] + [[['L'],['L']],[['nu'],['L']]])"
TChiChipmSlepL.conditionDescription = None
TChiChipmSlepL.condition = "cSim([[['L'],['L']],[['L'],['nu']]],[[['L'],['L']],[['nu'],['L']]]); cGtr([[['L'],['L']],[['nu'],['L']]],3.*[[['ta'],['ta']],[['nu'],['L']]]); cGtr([[['L'],['L']],[['L'],['nu']]],3.*[[['ta'],['ta']],[['L'],['nu']]]); cGtr([[['L'],['L']],[['nu'],['L']]],3.*[[['L'],['L']],[['nu'],['ta']]]); cGtr([[['L'],['L']],[['L'],['nu']]],3.*[[['L'],['L']],[['ta'],['nu']]]);cGtr([[['L'],['L']],[['nu'],['L']]],3.*[[['e'],['e']],[['nu'],['L']]]); cGtr([[['L'],['L']],[['L'],['nu']]],3.*[[['e'],['e']],[['L'],['nu']]]); cGtr([[['L'],['L']],[['nu'],['L']]],3.*[[['L'],['L']],[['nu'],['e']]]); cGtr([[['L'],['L']],[['L'],['nu']]],3.*[[['L'],['L']],[['e'],['nu']]])"
TChiChipmSlepL.source = "CMS"
#+++++++ next mass plane block ++++++++++++++
TChiChipmSlepL_1 = TChiChipmSlepL.addMassPlane(2*[[x, 0.05*x+0.95*y, y]])
TChiChipmSlepL_1.figure = 'Fig. 13'
TChiChipmSlepL_1.figureUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/Fig13_exclusion_TChiSlepSnu_2i_0_05.png'
TChiChipmSlepL_1.dataUrl = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13006'
TChiChipmSlepL_1.histoDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/exclusion_TChiSlepSnu_2i_0_05.root'
TChiChipmSlepL_1.dataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/exclusion_TChiSlepSnu_2i_0_05.root'
TChiChipmSlepL_1.exclusionDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/exclusion_TChiSlepSnu_2i_0_05.root'
TChiChipmSlepL_1.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/exclusion_TChiSlepSnu_2i_0_05.root', 'orig/exclusion_TChiSlepSnu_2i_0_05.root', 'orig/exclusion_TChiSlepSnu_2i_0_05.root', 'orig/exclusion_TChiSlepSnu_2i_0_05.root', 'orig/exclusion_TChiSlepSnu_2i_0_05.root', 'orig/exclusion_TChiSlepSnu_2i_0_05.root', 'orig/exclusion_TChiSlepSnu_2i_0_05.root'],
                 dataFormats= ['canvas', 'canvas', 'canvas', 'canvas', 'canvas', 'canvas', 'canvas'],objectNames= ['interpret', 'interpret', 'interpret', 'interpret', 'interpret', 'interpret', 'interpret'],indices= [9, 5, 6, 8, 8, 7, 2],units= [None, None, None, None, None, None, 'fb'])
#+++++++ next mass plane block ++++++++++++++
TChiChipmSlepL_2 = TChiChipmSlepL.addMassPlane(2*[[x, 0.5*x+0.5*y, y]])
TChiChipmSlepL_2.figure ='Fig. 12'
TChiChipmSlepL_2.figureUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/Fig13_exclusion_TChiSlepSnu_2i_0_5.png'
TChiChipmSlepL_2.dataUrl ='https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13006'
TChiChipmSlepL_2.histoDataUrl = 'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/exclusion_TChiSlepSnu_2i_0_5.root'
TChiChipmSlepL_2.dataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/exclusion_TChiSlepSnu_2i_0_5.root'
TChiChipmSlepL_2.exclusionDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/exclusion_TChiSlepSnu_2i_0_5.root'
TChiChipmSlepL_2.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/exclusion_TChiSlepSnu_2i_0_5.root', 'orig/exclusion_TChiSlepSnu_2i_0_5.root', 'orig/exclusion_TChiSlepSnu_2i_0_5.root', 'orig/exclusion_TChiSlepSnu_2i_0_5.root', 'orig/exclusion_TChiSlepSnu_2i_0_5.root', 'orig/exclusion_TChiSlepSnu_2i_0_5.root', 'orig/exclusion_TChiSlepSnu_2i_0_5.root'],
                 dataFormats= ['canvas', 'canvas', 'canvas', 'canvas', 'canvas', 'canvas', 'canvas'],objectNames= ['interpret', 'interpret', 'interpret', 'interpret', 'interpret', 'interpret', 'interpret'],indices= [8, 4, 5, 9, 7, 6, 2],units= [None, None, None, None, None, None, 'fb'])
#+++++++ next mass plane block ++++++++++++++
TChiChipmSlepL_3 = TChiChipmSlepL.addMassPlane(2*[[x, 0.95*x+0.05*y, y]])
TChiChipmSlepL_3.figure ='Fig. 14'
TChiChipmSlepL_3.figureUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/Fig13_exclusion_TChiSlepSnu_2i_0_95.png'
TChiChipmSlepL_3.dataUrl ='https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13006'
TChiChipmSlepL_3.histoDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/exclusion_TChiSlepSnu_2i_0_95.root'
TChiChipmSlepL_3.dataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/exclusion_TChiSlepSnu_2i_0_95.root'
TChiChipmSlepL_3.exclusionDataUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/exclusion_TChiSlepSnu_2i_0_95.root'
TChiChipmSlepL_3.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/exclusion_TChiSlepSnu_2i_0_95.root', 'orig/exclusion_TChiSlepSnu_2i_0_95.root', 'orig/exclusion_TChiSlepSnu_2i_0_95.root', 'orig/exclusion_TChiSlepSnu_2i_0_95.root', 'orig/exclusion_TChiSlepSnu_2i_0_95.root', 'orig/exclusion_TChiSlepSnu_2i_0_95.root', 'orig/exclusion_TChiSlepSnu_2i_0_95.root'],
                 dataFormats= ['canvas', 'canvas', 'canvas', 'canvas', 'canvas', 'canvas', 'canvas'],objectNames= ['interpret', 'interpret', 'interpret', 'interpret', 'interpret', 'interpret', 'interpret'],indices= [8, 4, 5, 11, 7, 6, 2],units= [None, None, None, None, None, None, 'fb'])



databaseCreator.create()
