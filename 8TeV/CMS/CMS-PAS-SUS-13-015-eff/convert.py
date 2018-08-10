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
info = MetaInfoInput('CMS-PAS-SUS-13-015')
info.url = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13015'
info.sqrts = 8
info.lumi = 19.4
info.prettyName = '>= 5(1b-)jets + Etmiss'
info.private = False
info.publication ='http://cds.cern.ch/record/1635353/files/SUS-13-015-pas.pdf'
info.comment ='PAS so no arxiv publication - only cds.cern; whenever the error on the expected SM background is asymmetric for some signal regions, the largest value is used'




#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("pTmiss200_Nb2")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "pTmiss200_Nb2", observedN = 83, expectedBG = 88.4 , bgError = 19.8, upperLimit = '2.068E+00*fb', expectedUpperLimit = '2.244E+00*fb')
#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.constraint = "[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition = None
T2tt.massConstraint = None
T2tt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.constraint ="[[['W','b']],[['W','b']]]"
T2ttoff.conditionDescription =None
T2ttoff.condition = None
T2ttoff.massConstraint = [['dm <= 169.']]*2
T2ttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.figure = 'Fig. 11'
T2tt_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13015/picked_eff_T2tt_baseline_2b.png"
T2tt_1.addSource('obsExclusion','orig/CMS-SUS-13-015_from_cMacro.root', 'root',   objectName = 'combined_obsExclOneTimesProspino')
T2tt_1.addSource('obsExclusionM1','orig/CMS-SUS-13-015_from_cMacro.root', 'root', objectName = 'combined_obsExclMinusSysErrProspino')
T2tt_1.addSource('obsExclusionP1','orig/CMS-SUS-13-015_from_cMacro.root', 'root', objectName = 'combined_obsExclPlusSysErrProspino')
T2tt_1.addSource('expExclusion','orig/CMS-SUS-13-015_from_cMacro.root', 'root',   objectName = 'combined_expExclOneTimesProspino')
T2tt_1.addSource('expExclusionM1','orig/CMS-SUS-13-015_from_cMacro.root', 'root', objectName = 'combined_obsExclMinusSysErrProspino')
T2tt_1.addSource('expExclusionP1','orig/CMS-SUS-13-015_from_cMacro.root', 'root', objectName = 'combined_obsExclPlusSysErrProspino')
T2tt_1.addSource('efficiencyMap','orig/eff_highMET.root', 'root', objectName = 'eff_highMET')
T2tt_1.dataUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13015"
T2ttoff.addMassPlane(T2tt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("pTmiss350_Nb1")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "pTmiss350_Nb1", observedN = 45, expectedBG = 40.9 , bgError = 9.6, upperLimit = '1.383E+00*fb', expectedUpperLimit = '1.189E+00*fb')
#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.constraint = "[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition = None
T2tt.massConstraint = None
T2tt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.constraint ="[[['W','b']],[['W','b']]]"
T2ttoff.conditionDescription =None
T2ttoff.condition = None
T2ttoff.massConstraint = [['dm <= 169.']]*2
T2ttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.figure = 'Fig. 11'
T2tt_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13015/picked_eff_T2tt_highMET.png"
T2tt_1.addSource('obsExclusion','orig/CMS-SUS-13-015_from_cMacro.root', 'root',   objectName = 'combined_obsExclOneTimesProspino')
T2tt_1.addSource('obsExclusionM1','orig/CMS-SUS-13-015_from_cMacro.root', 'root', objectName = 'combined_obsExclMinusSysErrProspino')
T2tt_1.addSource('obsExclusionP1','orig/CMS-SUS-13-015_from_cMacro.root', 'root', objectName = 'combined_obsExclPlusSysErrProspino')
T2tt_1.addSource('expExclusion','orig/CMS-SUS-13-015_from_cMacro.root', 'root',   objectName = 'combined_expExclOneTimesProspino')
T2tt_1.addSource('expExclusionM1','orig/CMS-SUS-13-015_from_cMacro.root', 'root', objectName = 'combined_obsExclMinusSysErrProspino')
T2tt_1.addSource('expExclusionP1','orig/CMS-SUS-13-015_from_cMacro.root', 'root', objectName = 'combined_obsExclPlusSysErrProspino')
T2tt_1.addSource('efficiencyMap','orig/eff_baseline_2b.root', 'root', objectName = 'eff_baseline_2b')
T2tt_1.dataUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13015"
T2ttoff.addMassPlane(T2tt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("pTmiss350_Nb2")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "pTmiss350_Nb2", observedN = 15, expectedBG = 8.6 , bgError = 7.1, upperLimit = '9.121E-01*fb', expectedUpperLimit = '5.572E-01*fb')
#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.constraint = "[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition = None
T2tt.massConstraint = None
T2tt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.constraint ="[[['W','b']],[['W','b']]]"
T2ttoff.conditionDescription =None
T2ttoff.condition = None
T2ttoff.massConstraint = [['dm <= 169.']]*2
T2ttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.figure = 'Fig. 11'
T2tt_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13015/picked_eff_T2tt_highMET_2b.png"
T2tt_1.addSource('obsExclusion','orig/CMS-SUS-13-015_from_cMacro.root', 'root',   objectName = 'combined_obsExclOneTimesProspino')
T2tt_1.addSource('obsExclusionM1','orig/CMS-SUS-13-015_from_cMacro.root', 'root', objectName = 'combined_obsExclMinusSysErrProspino')
T2tt_1.addSource('obsExclusionP1','orig/CMS-SUS-13-015_from_cMacro.root', 'root', objectName = 'combined_obsExclPlusSysErrProspino')
T2tt_1.addSource('expExclusion','orig/CMS-SUS-13-015_from_cMacro.root', 'root',   objectName = 'combined_expExclOneTimesProspino')
T2tt_1.addSource('expExclusionM1','orig/CMS-SUS-13-015_from_cMacro.root', 'root', objectName = 'combined_obsExclMinusSysErrProspino')
T2tt_1.addSource('expExclusionP1','orig/CMS-SUS-13-015_from_cMacro.root', 'root', objectName = 'combined_obsExclPlusSysErrProspino')
T2tt_1.addSource('efficiencyMap','orig/eff_highMET_2b.root', 'root', objectName = 'eff_highMET_2b')
T2tt_1.dataUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13015"
T2ttoff.addMassPlane(T2tt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("pTmiss200_Nb1")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "pTmiss200_Nb1", observedN = 254, expectedBG = 254.3 , bgError = 35, upperLimit = '3.892E+00*fb', expectedUpperLimit = '3.892E+00*fb')
#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.constraint = "[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition = None
T2tt.massConstraint = None
T2tt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.constraint ="[[['W','b']],[['W','b']]]"
T2ttoff.conditionDescription =None
T2ttoff.condition = None
T2ttoff.massConstraint = [['dm <= 169.']]*2
T2ttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.figure = 'Fig. 11'
T2tt_1.figureUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13015/picked_eff_T2tt_baseline.png"
T2tt_1.addSource('obsExclusion','orig/CMS-SUS-13-015_from_cMacro.root', 'root',   objectName = 'combined_obsExclOneTimesProspino')
T2tt_1.addSource('obsExclusionM1','orig/CMS-SUS-13-015_from_cMacro.root', 'root', objectName = 'combined_obsExclMinusSysErrProspino')
T2tt_1.addSource('obsExclusionP1','orig/CMS-SUS-13-015_from_cMacro.root', 'root', objectName = 'combined_obsExclPlusSysErrProspino')
T2tt_1.addSource('expExclusion','orig/CMS-SUS-13-015_from_cMacro.root', 'root',   objectName = 'combined_expExclOneTimesProspino')
T2tt_1.addSource('expExclusionM1','orig/CMS-SUS-13-015_from_cMacro.root', 'root', objectName = 'combined_obsExclMinusSysErrProspino')
T2tt_1.addSource('expExclusionP1','orig/CMS-SUS-13-015_from_cMacro.root', 'root', objectName = 'combined_obsExclPlusSysErrProspino')
T2tt_1.addSource('efficiencyMap','orig/eff_baseline.root', 'root', objectName = 'eff_baseline')
T2tt_1.dataUrl  = "https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13015"
T2ttoff.addMassPlane(T2tt_1)

databaseCreator.create()
