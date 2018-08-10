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
info = MetaInfoInput('CMS-SUS-13-006')
info.url ='https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13006'
info.sqrts = 8
info.lumi = 19.5
info.prettyName = 'EW productions with decays to leptons, W, Z, and Higgs'
info.private = False
info.arxiv = 'http://arxiv.org/abs/1405.7570'
info.contact = ''
info.publication = 'http://link.springer.com/article/10.1140%2Fepjc%2Fs10052-014-3036-7'
info.comment ='Using single lepton analysis EM'
info.supersedes ='CMS-PAS-SUS-12-022'




#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("MET_150")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "MET_150", observedN = 3, expectedBG = 3.8 , bgError = 1.0, upperLimit = '2.752E-01*fb', expectedUpperLimit = '2.698E-01*fb')
#+++++++ next txName block ++++++++++++++
TChiWH = dataset.addTxName('TChiWH')
TChiWH.constraint = "[[['W']],[['higgs']]]"
TChiWH.conditionDescription =None
TChiWH.condition =None
TChiWH.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
TChiWH_1 = TChiWH.addMassPlane([[x,y]]*2)
TChiWH_1.figure = 'Fig. 16(right)'
TChiWH_1.figureUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/Fig16_exclusion_TChiWH.png'
TChiWH_1.addSource('obsExclusion',"orig/exclusion_TChiWH.root", "canvas", objectName = "interpret", index = 9 )
TChiWH_1.addSource('obsExclusionM1',"orig/exclusion_TChiWH.root", "canvas", objectName = "interpret", index = 7 )
TChiWH_1.addSource('obsExclusionP1',"orig/exclusion_TChiWH.root", "canvas", objectName = "interpret", index = 6 )
TChiWH_1.addSource('expExclusion',"orig/exclusion_TChiWH.root", "canvas", objectName = "interpret", index = 8 )
TChiWH_1.addSource('expExclusionM1',"orig/exclusion_TChiWH.root", "canvas", objectName = "interpret", index = 4 )
TChiWH_1.addSource('expExclusionP1',"orig/exclusion_TChiWH.root", "canvas", objectName = "interpret", index = 5 )
TChiWH_1.addSource('efficiencyMap','orig/singlelep_results.root', 'root', objectName = 'h_eff_met150', index = None, scale = 0.01 )
TChiWH_1.dataUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/ss_eff_map.root"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("MET_175")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "MET_175", observedN = 3, expectedBG = 2.3 , bgError = 0.6, upperLimit = '3.036E-01*fb', expectedUpperLimit = '2.469E-01*fb')
#+++++++ next txName block ++++++++++++++
TChiWH = dataset.addTxName('TChiWH')
TChiWH.constraint = "[[['W']],[['higgs']]]"
TChiWH.conditionDescription =None
TChiWH.condition =None
TChiWH.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
TChiWH_1 = TChiWH.addMassPlane([[x,y]]*2)
TChiWH_1.figure = 'Fig. 16(right)'
TChiWH_1.figureUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/Fig16_exclusion_TChiWH.png'
TChiWH_1.addSource('obsExclusion',"orig/exclusion_TChiWH.root", "canvas", objectName = "interpret", index = 9 )
TChiWH_1.addSource('obsExclusionM1',"orig/exclusion_TChiWH.root", "canvas", objectName = "interpret", index = 7 )
TChiWH_1.addSource('obsExclusionP1',"orig/exclusion_TChiWH.root", "canvas", objectName = "interpret", index = 6 )
TChiWH_1.addSource('expExclusion',"orig/exclusion_TChiWH.root", "canvas", objectName = "interpret", index = 8 )
TChiWH_1.addSource('expExclusionM1',"orig/exclusion_TChiWH.root", "canvas", objectName = "interpret", index = 4 )
TChiWH_1.addSource('expExclusionP1',"orig/exclusion_TChiWH.root", "canvas", objectName = "interpret", index = 5 )
TChiWH_1.addSource('efficiencyMap','orig/singlelep_results.root', 'root', objectName = 'h_eff_met175', index = None, scale = 0.01 )
TChiWH_1.dataUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/ss_eff_map.root"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("MET_100")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "MET_100", observedN = 7, expectedBG = 7.7 , bgError = 1.9, upperLimit = '3.867E-01*fb', expectedUpperLimit = '3.883E-01*fb')
#+++++++ next txName block ++++++++++++++
TChiWH = dataset.addTxName('TChiWH')
TChiWH.constraint = "[[['W']],[['higgs']]]"
TChiWH.conditionDescription =None
TChiWH.condition =None
TChiWH.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
TChiWH_1 = TChiWH.addMassPlane([[x,y]]*2)
TChiWH_1.figure = 'Fig. 16(right)'
TChiWH_1.figureUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/Fig16_exclusion_TChiWH.png'
TChiWH_1.addSource('obsExclusion',"orig/exclusion_TChiWH.root", "canvas", objectName = "interpret", index = 9 )
TChiWH_1.addSource('obsExclusionM1',"orig/exclusion_TChiWH.root", "canvas", objectName = "interpret", index = 7 )
TChiWH_1.addSource('obsExclusionP1',"orig/exclusion_TChiWH.root", "canvas", objectName = "interpret", index = 6 )
TChiWH_1.addSource('expExclusion',"orig/exclusion_TChiWH.root", "canvas", objectName = "interpret", index = 8 )
TChiWH_1.addSource('expExclusionM1',"orig/exclusion_TChiWH.root", "canvas", objectName = "interpret", index = 4 )
TChiWH_1.addSource('expExclusionP1',"orig/exclusion_TChiWH.root", "canvas", objectName = "interpret", index = 5 )
TChiWH_1.addSource('efficiencyMap','orig/singlelep_results.root', 'root', objectName = 'h_eff_met100', index = None, scale = 0.01 )
TChiWH_1.dataUrl  = "https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13006"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("MET_125")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "MET_125", observedN = 6, expectedBG = 5.4 , bgError = 1.3, upperLimit = '3.887E-01*fb', expectedUpperLimit = '3.359E-01*fb')
#+++++++ next txName block ++++++++++++++
TChiWH = dataset.addTxName('TChiWH')
TChiWH.constraint = "[[['W']],[['higgs']]]"
TChiWH.conditionDescription =None
TChiWH.condition =None
TChiWH.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
TChiWH_1 = TChiWH.addMassPlane([[x,y]]*2)
TChiWH_1.figure = 'Fig. 16(right)'
TChiWH_1.figureUrl ='https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/Fig16_exclusion_TChiWH.png'
TChiWH_1.addSource('obsExclusion',"orig/exclusion_TChiWH.root", "canvas", objectName = "interpret", index = 9 )
TChiWH_1.addSource('obsExclusionM1',"orig/exclusion_TChiWH.root", "canvas", objectName = "interpret", index = 7 )
TChiWH_1.addSource('obsExclusionP1',"orig/exclusion_TChiWH.root", "canvas", objectName = "interpret", index = 6 )
TChiWH_1.addSource('expExclusion',"orig/exclusion_TChiWH.root", "canvas", objectName = "interpret", index = 8 )
TChiWH_1.addSource('expExclusionM1',"orig/exclusion_TChiWH.root", "canvas", objectName = "interpret", index = 4 )
TChiWH_1.addSource('expExclusionP1',"orig/exclusion_TChiWH.root", "canvas", objectName = "interpret", index = 5 )
TChiWH_1.addSource('efficiencyMap','orig/singlelep_results.root', 'root', objectName = 'h_eff_met125', index = None, scale = 0.01 )
TChiWH_1.dataUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13006/ss_eff_map.root"

databaseCreator.create()
