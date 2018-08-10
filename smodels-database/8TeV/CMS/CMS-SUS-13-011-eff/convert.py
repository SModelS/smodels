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
info = MetaInfoInput('CMS-SUS-13-011')
info.url = 'https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsSUS13011'
info.sqrts = 8
info.lumi = 19.5
info.prettyName = '1 lepton + >= 4 (1b-)jets + Etmiss'
info.private = False
info.arxiv = 'arXiv:1308.1586v2'
info.publication = 'Eur. Phys. J. C 73 (2013) 2677'




#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("HM250")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "HM250", observedN = 3, expectedBG = 9.5 , bgError = 2.8, upperLimit = '2.4696E-01*fb', expectedUpperLimit = '4.7010E-01*fb')
#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.constraint =  "[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition = None
T2tt.massConstraint = None
T2tt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.constraint = "[[['W','b']],[['W','b']]]"
T2ttoff.conditionDescription ="None"
T2ttoff.condition = None
T2ttoff.massConstraint = [['dm <= 169.0']]*2
T2ttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.figure = 'topneutralino_cutbased_efficiencies'
T2tt_1.figureUrl =  'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/topneutralino_cutbased_efficiencies.pdf'
T2tt_1.addSource('obsExclusion', 'orig/SUS13011_T2tt_exclusion.txt', 'txt')
T2tt_1.addSource('efficiencyMap','orig/topneutralino_cutbased_efficiencies.root', 'root', objectName = 'efficiency_HM250')
T2tt_1.dataUrl =  'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/topneutralino_cutbased_efficiencies.pdf'
T2ttoff.addMassPlane(T2tt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("HM150")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "HM150", observedN = 23, expectedBG = 29 , bgError = 7, upperLimit = '7.3407E-01*fb', expectedUpperLimit = '9.4090E-01*fb')
#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.constraint =  "[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition = None
T2tt.massConstraint = None
T2tt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.constraint = "[[['W','b']],[['W','b']]]"
T2ttoff.conditionDescription ="None"
T2ttoff.condition = None
T2ttoff.massConstraint = [['dm <= 169.0']]*2
T2ttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.figure = 'topneutralino_cutbased_efficiencies'
T2tt_1.figureUrl =  'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/topneutralino_cutbased_efficiencies.pdf'
T2tt_1.addSource('obsExclusion', 'orig/SUS13011_T2tt_exclusion.txt', 'txt')
T2tt_1.addSource('efficiencyMap','orig/topneutralino_cutbased_efficiencies.root', 'root', objectName = 'efficiency_HM150')
T2tt_1.dataUrl =  'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/topneutralino_cutbased_efficiencies.pdf'
T2ttoff.addMassPlane(T2tt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("HM300")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "HM300", observedN = 2, expectedBG = 4.7 , bgError = 1.4, upperLimit = '2.2546E-01*fb', expectedUpperLimit = '3.1259E-01*fb')
#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.constraint =  "[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition = None
T2tt.massConstraint = None
T2tt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.constraint = "[[['W','b']],[['W','b']]]"
T2ttoff.conditionDescription ="None"
T2ttoff.condition = None
T2ttoff.massConstraint = [['dm <= 169.0']]*2
T2ttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.figure = 'topneutralino_cutbased_efficiencies'
T2tt_1.figureUrl =  'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/topneutralino_cutbased_efficiencies.pdf'
T2tt_1.addSource('obsExclusion', 'orig/SUS13011_T2tt_exclusion.txt', 'txt')
T2tt_1.addSource('efficiencyMap','orig/topneutralino_cutbased_efficiencies.root', 'root', objectName = 'efficiency_HM300')
T2tt_1.dataUrl =  'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/topneutralino_cutbased_efficiencies.pdf'
T2ttoff.addMassPlane(T2tt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("LM300")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "LM300", observedN = 9, expectedBG = 11.5 , bgError = 3.6, upperLimit = '4.6036E-01*fb', expectedUpperLimit = '5.4228E-01*fb')
#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.constraint =  "[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition = None
T2tt.massConstraint = None
T2tt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.constraint = "[[['W','b']],[['W','b']]]"
T2ttoff.conditionDescription ="None"
T2ttoff.condition = None
T2ttoff.massConstraint = [['dm <= 169.0']]*2
T2ttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.figure = 'topneutralino_cutbased_efficiencies'
T2tt_1.figureUrl =  'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/topneutralino_cutbased_efficiencies.pdf'
T2tt_1.addSource('obsExclusion', 'orig/SUS13011_T2tt_exclusion.txt', 'txt')
T2tt_1.addSource('efficiencyMap','orig/topneutralino_cutbased_efficiencies.root', 'root', objectName = 'efficiency_LM300')
T2tt_1.dataUrl =  'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/topneutralino_cutbased_efficiencies.pdf'
T2ttoff.addMassPlane(T2tt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("HM200")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "HM200", observedN = 11, expectedBG = 17 , bgError = 5, upperLimit = '4.9064E-01*fb', expectedUpperLimit = '7.1149E-01*fb')
#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.constraint =  "[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition = None
T2tt.massConstraint = None
T2tt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.constraint = "[[['W','b']],[['W','b']]]"
T2ttoff.conditionDescription ="None"
T2ttoff.condition = None
T2ttoff.massConstraint = [['dm <= 169.0']]*2
T2ttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.figure = 'topneutralino_cutbased_efficiencies'
T2tt_1.figureUrl =  'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/topneutralino_cutbased_efficiencies.pdf'
T2tt_1.addSource('obsExclusion', 'orig/SUS13011_T2tt_exclusion.txt', 'txt')
T2tt_1.addSource('efficiencyMap','orig/topneutralino_cutbased_efficiencies.root', 'root', objectName = 'efficiency_HM200')
T2tt_1.dataUrl =  'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/topneutralino_cutbased_efficiencies.pdf'
T2ttoff.addMassPlane(T2tt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("LM200")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "LM200", observedN = 69, expectedBG = 83 , bgError = 21, upperLimit = '1.8797E+00*fb', expectedUpperLimit = '2.3503E+00*fb')
#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.constraint =  "[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition = None
T2tt.massConstraint = None
T2tt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.constraint = "[[['W','b']],[['W','b']]]"
T2ttoff.conditionDescription ="None"
T2ttoff.condition = None
T2ttoff.massConstraint = [['dm <= 169.0']]*2
T2ttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.figure = 'topneutralino_cutbased_efficiencies'
T2tt_1.figureUrl =  'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/topneutralino_cutbased_efficiencies.pdf'
T2tt_1.addSource('obsExclusion', 'orig/SUS13011_T2tt_exclusion.txt', 'txt')
T2tt_1.addSource('efficiencyMap','orig/topneutralino_cutbased_efficiencies.root', 'root', objectName = 'efficiency_LM200')
T2tt_1.dataUrl =  'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/topneutralino_cutbased_efficiencies.pdf'
T2ttoff.addMassPlane(T2tt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("LM150")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "LM150", observedN = 227, expectedBG = 251 , bgError = 50, upperLimit = '4.5368E+00*fb', expectedUpperLimit = '5.2995E+00*fb')
#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.constraint =  "[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition = None
T2tt.massConstraint = None
T2tt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.constraint = "[[['W','b']],[['W','b']]]"
T2ttoff.conditionDescription ="None"
T2ttoff.condition = None
T2ttoff.massConstraint = [['dm <= 169.0']]*2
T2ttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.figure = 'topneutralino_cutbased_efficiencies'
T2tt_1.figureUrl =  'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/topneutralino_cutbased_efficiencies.pdf'
T2tt_1.addSource('obsExclusion', 'orig/SUS13011_T2tt_exclusion.txt', 'txt')
T2tt_1.addSource('efficiencyMap','orig/topneutralino_cutbased_efficiencies.root', 'root', objectName = 'efficiency_LM150')
T2tt_1.dataUrl =  'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/topneutralino_cutbased_efficiencies.pdf'
T2ttoff.addMassPlane(T2tt_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("LM250")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "LM250", observedN = 21, expectedBG = 31 , bgError = 8, upperLimit = '7.0631E-01*fb', expectedUpperLimit = '1.0425E+00*fb')
#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.constraint =  "[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition = None
T2tt.massConstraint = None
T2tt.source = 'CMS'
#+++++++ next txName block ++++++++++++++
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.constraint = "[[['W','b']],[['W','b']]]"
T2ttoff.conditionDescription ="None"
T2ttoff.condition = None
T2ttoff.massConstraint = [['dm <= 169.0']]*2
T2ttoff.source = 'CMS'
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.figure = 'topneutralino_cutbased_efficiencies'
T2tt_1.figureUrl =  'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/topneutralino_cutbased_efficiencies.pdf'
T2tt_1.addSource('obsExclusion', 'orig/SUS13011_T2tt_exclusion.txt', 'txt')
T2tt_1.addSource('efficiencyMap','orig/topneutralino_cutbased_efficiencies.root', 'root', objectName = 'efficiency_LM250')
T2tt_1.dataUrl =  'https://twiki.cern.ch/twiki/pub/CMSPublic/PhysicsResultsSUS13011/topneutralino_cutbased_efficiencies.pdf'
T2ttoff.addMassPlane(T2tt_1)

databaseCreator.create()
