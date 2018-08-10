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
# databaseCreator.ncpus = 1

#+++++++ global info block ++++++++++++++
info = MetaInfoInput('ATLAS-SUSY-2013-04')
info.comment = 'Erratum: JHEP01(2014)109. T1tttt official efficiency maps from ATLAS collaboration;T1bbbb,T2tt,T1btbt,T5,T6bbWW,T5WW and T5ZZ efficiency maps created by the SModelS collaboration using MadAnalysis5'
info.sqrts = '8.0'
info.private = False
info.lumi = '20.3'
info.publication = 'http://link.springer.com/article/10.1007/JHEP10%282013%29130'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-04/'
# info.supersededBy =
info.arxiv = 'http://arxiv.org/abs/1308.1841'
info.contact = "ATLAS collaboration for T1tttt models;SModelS for T5WW and T5ZZ models"
info.prettyName = '0 leptons + >= 7-10 jets + Etmiss'
info.supersedes = 'ATLAS-CONF-2012-103'



#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("GtGrid_SR_8ij80_0bjet")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "GtGrid_SR_8ij80_0bjet", observedN = 2, expectedBG = 0.9 , bgError = 0.6)
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="None"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane([[x,y]]*2)
T1tttt.addSource('efficiencyMap',"./orig/GtGrid_SR_8ij80_0bjet.txt", "txt")
T1tttt.addSource('obsExclusion', "orig/exclusion_T1tttt.txt", "txt")
T1tttt.dataUrl="http://hepdata.cedar.ac.uk/resource/6095/bjetstream/GtGrid_SR_8ij80_0bjet.txt"
T1tttt.figure = 'Fig (aux). 11a'
T1tttt.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-04/figaux_11a.png'
#+++++++ next txName block ++++++++++++++
T5WW = dataset.addTxName('T5WW')
T5WW.constraint ="[[['jet','jet'],['W']],[['jet','jet'],['W']]]"
T5WW.conditionDescription =None
T5WW.condition =None
T5WW.massConstraint = None
T5WW.source = "SModelS"
#+++++++ next txName block ++++++++++++++
T5WWoff = dataset.addTxName('T5WWoff')
T5WWoff.constraint ="2.23*[[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
T5WWoff.conditionDescription =None
T5WWoff.condition =None
T5WWoff.massConstraint = [['dm >= 0.0','dm <= 76.']]*2
T5WWoff.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T5WW_1 = T5WW.addMassPlane([[x, 0.5*(x+y), y]]*2)
T5WW_1.addSource('efficiencyMap', 'orig/T5WW_x05/MA5_EM_T5WW_x05_>=8j80,0bjet.dat', 'txt')
T5WW_1.addSource('obsExclusion', 'orig/Exclusion_T5WW_x05.txt', 'txt')
T5WW_1.addSource('obsExclusionM1', 'orig/T5WW_Plus.txt', 'txt')
T5WW_1.addSource('obsExclusionP1', 'orig/T5WW_Minus.txt', 'txt')
T5WW_1.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5WW_2 = T5WW.addMassPlane([[x, 0.05*x + 0.95*y, y]]*2)
T5WW_2.addSource('efficiencyMap','orig/T5WW_x005/MA5_EM_T5WW_Glu005Neu095_>=8j80,0bjet.dat','txt')
T5WW_2.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5WW_3 = T5WW.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5WW_3.addSource('efficiencyMap','orig/T5WW_x095/MA5_EM_T5WW_Glu095Neu005_>=8j80,0bjet.dat','txt')
T5WW_3.dataUrl = None
T5WWoff.addMassPlane(T5WW_1)
T5WWoff.addMassPlane(T5WW_2)
T5WWoff.addMassPlane(T5WW_3)
#+++++++ next txName block ++++++++++++++
T5ZZ = dataset.addTxName('T5ZZ')
T5ZZ.constraint ="[[['jet','jet'],['Z']],[['jet','jet'],['Z']]]"
T5ZZ.conditionDescription =None
T5ZZ.condition =None
T5ZZ.massConstraint = None
T5ZZ.source = "SModelS"
#+++++++ next txName block ++++++++++++++
#+++++++ next mass plane block ++++++++++++++
T5ZZ_1 = T5ZZ.addMassPlane([[x,0.5*(x+y),y]]*2)
T5ZZ_1.addSource('efficiencyMap','orig/T5ZZ_x05/MA5_EM_T5ZZ_x05_>=8j80,0bjet.dat', 'txt')
T5ZZ_1.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5ZZ_2 = T5ZZ.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
T5ZZ_2.addSource('efficiencyMap','orig/T5ZZ_x005/MA5_EM_T5ZZ_Glu005Neu095_>=8j80,0bjet.dat', 'txt')
T5ZZ_2.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5ZZ_3 = T5ZZ.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5ZZ_3.addSource('efficiencyMap','orig/T5ZZ_x095/MA5_EM_T5ZZ_Glu095Neu005_>=8j80,0bjet.dat', 'txt')
T5ZZ_3.dataUrl = None
"""
T5ZZoff = dataset.addTxName('T5ZZoff')
T5ZZoff.constraint ="2.1*[[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
T5ZZoff.conditionDescription =None
T5ZZoff.condition =None
T5ZZoff.massConstraint = [['dm >= 0.0','dm <= 86.']]*2
T5ZZoff.source = "SModelS"
T5ZZoff.addMassPlane(T5ZZ_1)
T5ZZoff.addMassPlane(T5ZZ_2)
T5ZZoff.addMassPlane(T5ZZ_3)
"""
#+++++++ next mass plane block ++++++++++++++
fullpath = "./orig/atlas_susy_2013_04_T2tt_1_EM_MAPS/MA5_EM_T2tt_1_>=8j80,0bjet.dat"
T2tt = dataset.addTxName('T2tt')
T2tt.checked ="None"
T2tt.constraint ="[[['t']],[['t']]]"
T2tt.massConstraint = None
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = "SModelS"
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.addSource('efficiencyMap', fullpath , "txt")
T2tt_1.dataUrl=None
#+++++++ next mass plane block ++++++++++++++
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.checked ="None"
T2ttoff.constraint ="[[['W','b']],[['W','b']]]"
T2ttoff.conditionDescription ="None"
T2ttoff.condition ="None"
T2ttoff.source = "SModelS"
T2ttoff_1 = T2ttoff.addMassPlane([[x,y]]*2)
T2ttoff.massConstraint = [['dm <= 169.0']]*2
T2ttoff_1.addSource('efficiencyMap', fullpath , "txt")
T2ttoff_1.dataUrl=None
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="None"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb.addSource('efficiencyMap',"./orig/atlas_susy_2013_04_T1bbbb_1_EM_MAPS/MA5_EM_T1bbbb_1_>=8j80,0bjet.dat", "txt")
T1bbbb.dataUrl= None
T1btbt = dataset.addTxName('T1btbt')
T1btbt.checked ="None"
T1btbt.constraint ="[[['b','t']],[['b','t']]]"
T1btbt.conditionDescription ="None"
T1btbt.condition ="None"
T1btbt.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T1btbt = T1btbt.addMassPlane([[x,y]]*2)
T1btbt.addSource('efficiencyMap',"./orig/atlas_susy_2013_04_T1btbt_1_EM_MAPS/MA5_EM_T1btbt_1_>=8j80,0bjet.dat", "txt")
T1btbt.dataUrl= None
T5 = dataset.addTxName('T5')
T5.constraint ="[[['jet'],['jet']],[['jet'],['jet']]]"
T5.conditionDescription =None
T5.condition =None
T5.massConstraint = None
T5.source = "SModelS"
## T5.massConstraint = [['dm >= 0.0','dm > 76.']]*2
T5.dataUrl = None
#+++++++ next txName block ++++++++++++++
#+++++++ next mass plane block ++++++++++++++
T5_1 = T5.addMassPlane([[x,0.5*(x+y),y]]*2)
T5_1.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x05_EM_MAPS/MA5_EM_T5_x05_GtGrid_SR_8ij80_0bjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T5_2 = T5.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
T5_2.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x005_EM_MAPS/MA5_EM_T5_x005_GtGrid_SR_8ij80_0bjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T5_3 = T5.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5_3.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x095_EM_MAPS/MA5_EM_T5_x095_GtGrid_SR_8ij80_0bjet.dat', 'txt')
T6bbWW = dataset.addTxName('T6bbWW')
T6bbWW.constraint = "[[['b'],['W']],[['b'],['W']]]"
T6bbWW.conditionDescription =None
T6bbWW.condition =None
T6bbWW.massConstraint = None
T6bbWW.source = "SModelS"
T6bbWW.dataUrl = None
T6bbWW_1 = T6bbWW.addMassPlane([[x, 0.1*x + 0.9*y, y]]*2)
T6bbWW_1.addSource('efficiencyMap', 'orig/atlas_susy_2013_04_T6bbWW_x01_EM_MAPS/MA5_EM_T6bbWW_x01_GtGrid_SR_8ij80_0bjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T6bbWW_2 = T6bbWW.addMassPlane([[x, 0.5*x + 0.5*y, y]]*2)
T6bbWW_2.addSource('efficiencyMap','orig/atlas_susy_2013_04_T6bbWW_x05_EM_MAPS/MA5_EM_T6bbWW_x05_GtGrid_SR_8ij80_0bjet.dat','txt')
#+++++++ next mass plane block ++++++++++++++
T6bbWW_3 = T6bbWW.addMassPlane([[x, 0.9*x + 0.1*y, y]]*2)
T6bbWW_3.addSource('efficiencyMap','orig/atlas_susy_2013_04_T6bbWW_x09_EM_MAPS/MA5_EM_T6bbWW_x09_GtGrid_SR_8ij80_0bjet.dat','txt')


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("GtGrid_SR_9ej50_0bjet")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "GtGrid_SR_9ej50_0bjet", observedN = 5, expectedBG = 3.3 , bgError = 0.7)
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="None"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane([[x,y]]*2)
T1tttt.addSource('efficiencyMap',"./orig/GtGrid_SR_9ej50_0bjet.txt", "txt")
T1tttt.addSource('obsExclusion', "orig/exclusion_T1tttt.txt", "txt")
T1tttt.dataUrl="http://hepdata.cedar.ac.uk/resource/6095/bjetstream/GtGrid_SR_9ej50_0bjet.txt"
T1tttt.figure = 'Fig (aux). 11a'
T1tttt.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-04/figaux_11a.png'
#+++++++ next txName block ++++++++++++++
T5WW = dataset.addTxName('T5WW')
T5WW.constraint ="[[['jet','jet'],['W']],[['jet','jet'],['W']]]"
T5WW.conditionDescription =None
T5WW.condition =None
T5WW.massConstraint = None
T5WW.source = "SModelS"
#+++++++ next txName block ++++++++++++++
T5WWoff = dataset.addTxName('T5WWoff')
T5WWoff.constraint ="2.23*[[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
T5WWoff.conditionDescription =None
T5WWoff.condition =None
T5WWoff.massConstraint = [['dm >= 0.0','dm <= 76.']]*2
T5WWoff.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T5WW_1 = T5WW.addMassPlane([[x, 0.5*(x+y), y]]*2)
T5WW_1.addSource('efficiencyMap', 'orig/T5WW_x05/MA5_EM_T5WW_x05_9j50,0bjet.dat', 'txt')
T5WW_1.addSource('obsExclusion', 'orig/Exclusion_T5WW_x05.txt', 'txt')
T5WW_1.addSource('obsExclusionM1', 'orig/T5WW_Plus.txt', 'txt')
T5WW_1.addSource('obsExclusionP1', 'orig/T5WW_Minus.txt', 'txt')
T5WW_1.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5WW_2 = T5WW.addMassPlane([[x, 0.05*x + 0.95*y, y]]*2)
T5WW_2.addSource('efficiencyMap','orig/T5WW_x005/MA5_EM_T5WW_Glu005Neu095_9j50,0bjet.dat','txt')
T5WW_2.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5WW_3 = T5WW.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5WW_3.addSource('efficiencyMap','orig/T5WW_x095/MA5_EM_T5WW_Glu095Neu005_9j50,0bjet.dat','txt')
T5WW_3.dataUrl = None
T5WWoff.addMassPlane(T5WW_1)
T5WWoff.addMassPlane(T5WW_2)
T5WWoff.addMassPlane(T5WW_3)
#+++++++ next txName block ++++++++++++++
T5ZZ = dataset.addTxName('T5ZZ')
T5ZZ.constraint ="[[['jet','jet'],['Z']],[['jet','jet'],['Z']]]"
T5ZZ.conditionDescription =None
T5ZZ.condition =None
T5ZZ.massConstraint = None
T5ZZ.source = "SModelS"
#+++++++ next txName block ++++++++++++++
#+++++++ next mass plane block ++++++++++++++
T5ZZ_1 = T5ZZ.addMassPlane([[x,0.5*(x+y),y]]*2)
T5ZZ_1.addSource('efficiencyMap','orig/T5ZZ_x05/MA5_EM_T5ZZ_x05_9j50,0bjet.dat', 'txt')
T5ZZ_1.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5ZZ_2 = T5ZZ.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
T5ZZ_2.addSource('efficiencyMap','orig/T5ZZ_x005/MA5_EM_T5ZZ_Glu005Neu095_9j50,0bjet.dat', 'txt')
T5ZZ_2.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5ZZ_3 = T5ZZ.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5ZZ_3.addSource('efficiencyMap','orig/T5ZZ_x095/MA5_EM_T5ZZ_Glu095Neu005_9j50,0bjet.dat', 'txt')
T5ZZ_3.dataUrl = None
"""
T5ZZoff = dataset.addTxName('T5ZZoff')
T5ZZoff.constraint ="2.1*[[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
T5ZZoff.conditionDescription =None
T5ZZoff.condition =None
T5ZZoff.massConstraint = [['dm >= 0.0','dm <= 86.']]*2
T5ZZoff.source = "SModelS"
T5ZZoff.addMassPlane(T5ZZ_1)
T5ZZoff.addMassPlane(T5ZZ_2)
T5ZZoff.addMassPlane(T5ZZ_3)
"""
#+++++++ next mass plane block ++++++++++++++
fullpath = "./orig/atlas_susy_2013_04_T2tt_1_EM_MAPS/MA5_EM_T2tt_1_9j50,0bjet.dat"
T2tt = dataset.addTxName('T2tt')
T2tt.checked ="None"
T2tt.constraint ="[[['t']],[['t']]]"
T2tt.massConstraint = None
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = "SModelS"
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.addSource('efficiencyMap', fullpath , "txt")
T2tt_1.dataUrl=None
#+++++++ next mass plane block ++++++++++++++
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.checked ="None"
T2ttoff.constraint ="[[['W','b']],[['W','b']]]"
T2ttoff.conditionDescription ="None"
T2ttoff.condition ="None"
T2ttoff.source = "SModelS"
T2ttoff_1 = T2ttoff.addMassPlane([[x,y]]*2)
T2ttoff.massConstraint = [['dm <= 169.0']]*2
T2ttoff_1.addSource('efficiencyMap', fullpath , "txt")
T2ttoff_1.dataUrl=None
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="None"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb.addSource('efficiencyMap',"./orig/atlas_susy_2013_04_T1bbbb_1_EM_MAPS/MA5_EM_T1bbbb_1_9j50,0bjet.dat", "txt")
T1bbbb.dataUrl= None
T1btbt = dataset.addTxName('T1btbt')
T1btbt.checked ="None"
T1btbt.constraint ="[[['b','t']],[['b','t']]]"
T1btbt.conditionDescription ="None"
T1btbt.condition ="None"
T1btbt.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T1btbt = T1btbt.addMassPlane([[x,y]]*2)
T1btbt.addSource('efficiencyMap',"./orig/atlas_susy_2013_04_T1btbt_1_EM_MAPS/MA5_EM_T1btbt_1_9j50,0bjet.dat", "txt")
T1btbt.dataUrl= None
T5 = dataset.addTxName('T5')
T5.constraint ="[[['jet'],['jet']],[['jet'],['jet']]]"
T5.conditionDescription =None
T5.condition =None
T5.massConstraint = None
T5.source = "SModelS"
## T5.massConstraint = [['dm >= 0.0','dm > 76.']]*2
T5.dataUrl = None
#+++++++ next txName block ++++++++++++++
#+++++++ next mass plane block ++++++++++++++
T5_1 = T5.addMassPlane([[x,0.5*(x+y),y]]*2)
T5_1.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x05_EM_MAPS/MA5_EM_T5_x05_GtGrid_SR_9ej50_0bjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T5_2 = T5.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
T5_2.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x005_EM_MAPS/MA5_EM_T5_x005_GtGrid_SR_9ej50_0bjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T5_3 = T5.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5_3.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x095_EM_MAPS/MA5_EM_T5_x095_GtGrid_SR_9ej50_0bjet.dat', 'txt')
T6bbWW = dataset.addTxName('T6bbWW')
T6bbWW.constraint = "[[['b'],['W']],[['b'],['W']]]"
T6bbWW.conditionDescription =None
T6bbWW.condition =None
T6bbWW.massConstraint = None
T6bbWW.source = "SModelS"
T6bbWW.dataUrl = None
T6bbWW_1 = T6bbWW.addMassPlane([[x, 0.1*x + 0.9*y, y]]*2)
T6bbWW_1.addSource('efficiencyMap', 'orig/atlas_susy_2013_04_T6bbWW_x01_EM_MAPS/MA5_EM_T6bbWW_x01_GtGrid_SR_9ej50_0bjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T6bbWW_2 = T6bbWW.addMassPlane([[x, 0.5*x + 0.5*y, y]]*2)
T6bbWW_2.addSource('efficiencyMap','orig/atlas_susy_2013_04_T6bbWW_x05_EM_MAPS/MA5_EM_T6bbWW_x05_GtGrid_SR_9ej50_0bjet.dat','txt')
#+++++++ next mass plane block ++++++++++++++
T6bbWW_3 = T6bbWW.addMassPlane([[x, 0.9*x + 0.1*y, y]]*2)
T6bbWW_3.addSource('efficiencyMap','orig/atlas_susy_2013_04_T6bbWW_x09_EM_MAPS/MA5_EM_T6bbWW_x09_GtGrid_SR_9ej50_0bjet.dat','txt')


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("GtGrid_SR_8ej50_0bjet")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "GtGrid_SR_8ej50_0bjet", observedN = 40, expectedBG = 35 , bgError = 4)
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="None"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane([[x,y]]*2)
T1tttt.addSource('efficiencyMap',"./orig/GtGrid_SR_8ej50_0bjet.txt", "txt")
T1tttt.addSource('obsExclusion', "orig/exclusion_T1tttt.txt", "txt")
T1tttt.dataUrl="http://hepdata.cedar.ac.uk/resource/6095/bjetstream/GtGrid_SR_8ej50_0bjet.txt"
T1tttt.figure = 'Fig (aux). 11a'
T1tttt.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-04/figaux_11a.png'
#+++++++ next txName block ++++++++++++++
T5WW = dataset.addTxName('T5WW')
T5WW.constraint ="[[['jet','jet'],['W']],[['jet','jet'],['W']]]"
T5WW.conditionDescription =None
T5WW.condition =None
T5WW.massConstraint = None
T5WW.source = "SModelS"
#+++++++ next txName block ++++++++++++++
T5WWoff = dataset.addTxName('T5WWoff')
T5WWoff.constraint ="2.23*[[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
T5WWoff.conditionDescription =None
T5WWoff.condition =None
T5WWoff.massConstraint = [['dm >= 0.0','dm <= 76.']]*2
T5WWoff.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T5WW_1 = T5WW.addMassPlane([[x, 0.5*(x+y), y]]*2)
T5WW_1.addSource('efficiencyMap', 'orig/T5WW_x05/MA5_EM_T5WW_x05_8j50,0bjet.dat', 'txt')
T5WW_1.addSource('obsExclusion', 'orig/Exclusion_T5WW_x05.txt', 'txt')
T5WW_1.addSource('obsExclusionM1', 'orig/T5WW_Plus.txt', 'txt')
T5WW_1.addSource('obsExclusionP1', 'orig/T5WW_Minus.txt', 'txt')
T5WW_1.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5WW_2 = T5WW.addMassPlane([[x, 0.05*x + 0.95*y, y]]*2)
T5WW_2.addSource('efficiencyMap','orig/T5WW_x005/MA5_EM_T5WW_Glu005Neu095_8j50,0bjet.dat','txt')
T5WW_2.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5WW_3 = T5WW.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5WW_3.addSource('efficiencyMap','orig/T5WW_x095/MA5_EM_T5WW_Glu095Neu005_8j50,0bjet.dat','txt')
T5WW_3.dataUrl = None
T5WWoff.addMassPlane(T5WW_1)
T5WWoff.addMassPlane(T5WW_2)
T5WWoff.addMassPlane(T5WW_3)
#+++++++ next txName block ++++++++++++++
T5ZZ = dataset.addTxName('T5ZZ')
T5ZZ.constraint ="[[['jet','jet'],['Z']],[['jet','jet'],['Z']]]"
T5ZZ.conditionDescription =None
T5ZZ.condition =None
T5ZZ.massConstraint = None
T5ZZ.source = "SModelS"
#+++++++ next txName block ++++++++++++++
#+++++++ next mass plane block ++++++++++++++
T5ZZ_1 = T5ZZ.addMassPlane([[x,0.5*(x+y),y]]*2)
T5ZZ_1.addSource('efficiencyMap','orig/T5ZZ_x05/MA5_EM_T5ZZ_x05_8j50,0bjet.dat', 'txt')
T5ZZ_1.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5ZZ_2 = T5ZZ.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
T5ZZ_2.addSource('efficiencyMap','orig/T5ZZ_x005/MA5_EM_T5ZZ_Glu005Neu095_8j50,0bjet.dat', 'txt')
T5ZZ_2.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5ZZ_3 = T5ZZ.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5ZZ_3.addSource('efficiencyMap','orig/T5ZZ_x095/MA5_EM_T5ZZ_Glu095Neu005_8j50,0bjet.dat', 'txt')
T5ZZ_3.dataUrl = None
"""
T5ZZoff = dataset.addTxName('T5ZZoff')
T5ZZoff.constraint ="2.1*[[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
T5ZZoff.conditionDescription =None
T5ZZoff.condition =None
T5ZZoff.massConstraint = [['dm >= 0.0','dm <= 86.']]*2
T5ZZoff.source = "SModelS"
T5ZZoff.addMassPlane(T5ZZ_1)
T5ZZoff.addMassPlane(T5ZZ_2)
T5ZZoff.addMassPlane(T5ZZ_3)
"""
#+++++++ next mass plane block ++++++++++++++
fullpath = "./orig/atlas_susy_2013_04_T2tt_1_EM_MAPS/MA5_EM_T2tt_1_8j50,0bjet.dat"
T2tt = dataset.addTxName('T2tt')
T2tt.checked ="None"
T2tt.constraint ="[[['t']],[['t']]]"
T2tt.massConstraint = None
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = "SModelS"
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.addSource('efficiencyMap', fullpath , "txt")
T2tt_1.dataUrl=None
#+++++++ next mass plane block ++++++++++++++
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.checked ="None"
T2ttoff.constraint ="[[['W','b']],[['W','b']]]"
T2ttoff.conditionDescription ="None"
T2ttoff.condition ="None"
T2ttoff.source = "SModelS"
T2ttoff_1 = T2ttoff.addMassPlane([[x,y]]*2)
T2ttoff.massConstraint = [['dm <= 169.0']]*2
T2ttoff_1.addSource('efficiencyMap', fullpath , "txt")
T2ttoff_1.dataUrl=None
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="None"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb.addSource('efficiencyMap',"./orig/atlas_susy_2013_04_T1bbbb_1_EM_MAPS/MA5_EM_T1bbbb_1_8j50,0bjet.dat", "txt")
T1bbbb.dataUrl= None
T1btbt = dataset.addTxName('T1btbt')
T1btbt.checked ="None"
T1btbt.constraint ="[[['b','t']],[['b','t']]]"
T1btbt.conditionDescription ="None"
T1btbt.condition ="None"
T1btbt.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T1btbt = T1btbt.addMassPlane([[x,y]]*2)
T1btbt.addSource('efficiencyMap',"./orig/atlas_susy_2013_04_T1btbt_1_EM_MAPS/MA5_EM_T1btbt_1_8j50,0bjet.dat", "txt")
T1btbt.dataUrl= None
T5 = dataset.addTxName('T5')
T5.constraint ="[[['jet'],['jet']],[['jet'],['jet']]]"
T5.conditionDescription =None
T5.condition =None
T5.massConstraint = None
T5.source = "SModelS"
## T5.massConstraint = [['dm >= 0.0','dm > 76.']]*2
T5.dataUrl = None
#+++++++ next txName block ++++++++++++++
#+++++++ next mass plane block ++++++++++++++
T5_1 = T5.addMassPlane([[x,0.5*(x+y),y]]*2)
T5_1.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x05_EM_MAPS/MA5_EM_T5_x05_GtGrid_SR_8ej50_0bjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T5_2 = T5.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
T5_2.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x005_EM_MAPS/MA5_EM_T5_x005_GtGrid_SR_8ej50_0bjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T5_3 = T5.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5_3.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x095_EM_MAPS/MA5_EM_T5_x095_GtGrid_SR_8ej50_0bjet.dat', 'txt')
T6bbWW = dataset.addTxName('T6bbWW')
T6bbWW.constraint = "[[['b'],['W']],[['b'],['W']]]"
T6bbWW.conditionDescription =None
T6bbWW.condition =None
T6bbWW.massConstraint = None
T6bbWW.source = "SModelS"
T6bbWW.dataUrl = None
T6bbWW_1 = T6bbWW.addMassPlane([[x, 0.1*x + 0.9*y, y]]*2)
T6bbWW_1.addSource('efficiencyMap', 'orig/atlas_susy_2013_04_T6bbWW_x01_EM_MAPS/MA5_EM_T6bbWW_x01_GtGrid_SR_8ej50_0bjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T6bbWW_2 = T6bbWW.addMassPlane([[x, 0.5*x + 0.5*y, y]]*2)
T6bbWW_2.addSource('efficiencyMap','orig/atlas_susy_2013_04_T6bbWW_x05_EM_MAPS/MA5_EM_T6bbWW_x05_GtGrid_SR_8ej50_0bjet.dat','txt')
#+++++++ next mass plane block ++++++++++++++
T6bbWW_3 = T6bbWW.addMassPlane([[x, 0.9*x + 0.1*y, y]]*2)
T6bbWW_3.addSource('efficiencyMap','orig/atlas_susy_2013_04_T6bbWW_x09_EM_MAPS/MA5_EM_T6bbWW_x09_GtGrid_SR_8ej50_0bjet.dat','txt')


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("GtGrid_SR_8ij80_2ibjet")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "GtGrid_SR_8ij80_2ibjet", observedN = 3, expectedBG = 3.3 , bgError = 2.2)
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="None"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane([[x,y]]*2)
T1tttt.addSource('efficiencyMap',"./orig/GtGrid_SR_8ij80_2ibjet.txt", "txt")
T1tttt.addSource('obsExclusion', "orig/exclusion_T1tttt.txt", "txt")
T1tttt.dataUrl="http://hepdata.cedar.ac.uk/resource/6095/bjetstream/GtGrid_SR_8ij80_2ibjet.txt"
T1tttt.figure = 'Fig (aux). 11a'
T1tttt.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-04/figaux_11a.png'
#+++++++ next txName block ++++++++++++++
T5WW = dataset.addTxName('T5WW')
T5WW.constraint ="[[['jet','jet'],['W']],[['jet','jet'],['W']]]"
T5WW.conditionDescription =None
T5WW.condition =None
T5WW.massConstraint = None
T5WW.source = "SModelS"
#+++++++ next txName block ++++++++++++++
T5WWoff = dataset.addTxName('T5WWoff')
T5WWoff.constraint ="2.23*[[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
T5WWoff.conditionDescription =None
T5WWoff.condition =None
T5WWoff.massConstraint = [['dm >= 0.0','dm <= 76.']]*2
T5WWoff.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T5WW_1 = T5WW.addMassPlane([[x, 0.5*(x+y), y]]*2)
T5WW_1.addSource('efficiencyMap', 'orig/T5WW_x05/MA5_EM_T5WW_x05_>=8j80,>=2bjets.dat', 'txt')
T5WW_1.addSource('obsExclusion', 'orig/Exclusion_T5WW_x05.txt', 'txt')
T5WW_1.addSource('obsExclusionM1', 'orig/T5WW_Plus.txt', 'txt')
T5WW_1.addSource('obsExclusionP1', 'orig/T5WW_Minus.txt', 'txt')
T5WW_1.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5WW_2 = T5WW.addMassPlane([[x, 0.05*x + 0.95*y, y]]*2)
T5WW_2.addSource('efficiencyMap','orig/T5WW_x005/MA5_EM_T5WW_Glu005Neu095_>=8j80,>=2bjets.dat','txt')
T5WW_2.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5WW_3 = T5WW.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5WW_3.addSource('efficiencyMap','orig/T5WW_x095/MA5_EM_T5WW_Glu095Neu005_>=8j80,>=2bjets.dat','txt')
T5WW_3.dataUrl = None
T5WWoff.addMassPlane(T5WW_1)
T5WWoff.addMassPlane(T5WW_2)
T5WWoff.addMassPlane(T5WW_3)
#+++++++ next txName block ++++++++++++++
T5ZZ = dataset.addTxName('T5ZZ')
T5ZZ.constraint ="[[['jet','jet'],['Z']],[['jet','jet'],['Z']]]"
T5ZZ.conditionDescription =None
T5ZZ.condition =None
T5ZZ.massConstraint = None
T5ZZ.source = "SModelS"
#+++++++ next txName block ++++++++++++++
#+++++++ next mass plane block ++++++++++++++
T5ZZ_1 = T5ZZ.addMassPlane([[x,0.5*(x+y),y]]*2)
T5ZZ_1.addSource('efficiencyMap','orig/T5ZZ_x05/MA5_EM_T5ZZ_x05_>=8j80,>=2bjets.dat', 'txt')
T5ZZ_1.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5ZZ_2 = T5ZZ.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
T5ZZ_2.addSource('efficiencyMap','orig/T5ZZ_x005/MA5_EM_T5ZZ_Glu005Neu095_>=8j80,>=2bjets.dat', 'txt')
T5ZZ_2.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5ZZ_3 = T5ZZ.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5ZZ_3.addSource('efficiencyMap','orig/T5ZZ_x095/MA5_EM_T5ZZ_Glu095Neu005_>=8j80,>=2bjets.dat', 'txt')
T5ZZ_3.dataUrl = None
"""
T5ZZoff = dataset.addTxName('T5ZZoff')
T5ZZoff.constraint ="2.1*[[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
T5ZZoff.conditionDescription =None
T5ZZoff.condition =None
T5ZZoff.massConstraint = [['dm >= 0.0','dm <= 86.']]*2
T5ZZoff.source = "SModelS"
T5ZZoff.addMassPlane(T5ZZ_1)
T5ZZoff.addMassPlane(T5ZZ_2)
T5ZZoff.addMassPlane(T5ZZ_3)
"""
#+++++++ next mass plane block ++++++++++++++
fullpath = "./orig/atlas_susy_2013_04_T2tt_1_EM_MAPS/MA5_EM_T2tt_1_>=8j80,>=2bjets.dat"
T2tt = dataset.addTxName('T2tt')
T2tt.checked ="None"
T2tt.constraint ="[[['t']],[['t']]]"
T2tt.massConstraint = None
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = "SModelS"
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.addSource('efficiencyMap', fullpath , "txt")
T2tt_1.dataUrl=None
#+++++++ next mass plane block ++++++++++++++
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.checked ="None"
T2ttoff.constraint ="[[['W','b']],[['W','b']]]"
T2ttoff.conditionDescription ="None"
T2ttoff.condition ="None"
T2ttoff.source = "SModelS"
T2ttoff_1 = T2ttoff.addMassPlane([[x,y]]*2)
T2ttoff.massConstraint = [['dm <= 169.0']]*2
T2ttoff_1.addSource('efficiencyMap', fullpath , "txt")
T2ttoff_1.dataUrl=None
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="None"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb.addSource('efficiencyMap',"./orig/atlas_susy_2013_04_T1bbbb_1_EM_MAPS/MA5_EM_T1bbbb_1_>=8j80,>=2bjets.dat", "txt")
T1bbbb.dataUrl= None
T1btbt = dataset.addTxName('T1btbt')
T1btbt.checked ="None"
T1btbt.constraint ="[[['b','t']],[['b','t']]]"
T1btbt.conditionDescription ="None"
T1btbt.condition ="None"
T1btbt.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T1btbt = T1btbt.addMassPlane([[x,y]]*2)
T1btbt.addSource('efficiencyMap',"./orig/atlas_susy_2013_04_T1btbt_1_EM_MAPS/MA5_EM_T1btbt_1_>=8j80,>=2bjets.dat", "txt")
T1btbt.dataUrl= None
T5 = dataset.addTxName('T5')
T5.constraint ="[[['jet'],['jet']],[['jet'],['jet']]]"
T5.conditionDescription =None
T5.condition =None
T5.massConstraint = None
T5.source = "SModelS"
## T5.massConstraint = [['dm >= 0.0','dm > 76.']]*2
T5.dataUrl = None
#+++++++ next txName block ++++++++++++++
#+++++++ next mass plane block ++++++++++++++
T5_1 = T5.addMassPlane([[x,0.5*(x+y),y]]*2)
T5_1.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x05_EM_MAPS/MA5_EM_T5_x05_GtGrid_SR_8ij80_2ibjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T5_2 = T5.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
T5_2.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x005_EM_MAPS/MA5_EM_T5_x005_GtGrid_SR_8ij80_2ibjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T5_3 = T5.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5_3.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x095_EM_MAPS/MA5_EM_T5_x095_GtGrid_SR_8ij80_2ibjet.dat', 'txt')
T6bbWW = dataset.addTxName('T6bbWW')
T6bbWW.constraint = "[[['b'],['W']],[['b'],['W']]]"
T6bbWW.conditionDescription =None
T6bbWW.condition =None
T6bbWW.massConstraint = None
T6bbWW.source = "SModelS"
T6bbWW.dataUrl = None
T6bbWW_1 = T6bbWW.addMassPlane([[x, 0.1*x + 0.9*y, y]]*2)
T6bbWW_1.addSource('efficiencyMap', 'orig/atlas_susy_2013_04_T6bbWW_x01_EM_MAPS/MA5_EM_T6bbWW_x01_GtGrid_SR_8ij80_2ibjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T6bbWW_2 = T6bbWW.addMassPlane([[x, 0.5*x + 0.5*y, y]]*2)
T6bbWW_2.addSource('efficiencyMap','orig/atlas_susy_2013_04_T6bbWW_x05_EM_MAPS/MA5_EM_T6bbWW_x05_GtGrid_SR_8ij80_2ibjet.dat','txt')
#+++++++ next mass plane block ++++++++++++++
T6bbWW_3 = T6bbWW.addMassPlane([[x, 0.9*x + 0.1*y, y]]*2)
T6bbWW_3.addSource('efficiencyMap','orig/atlas_susy_2013_04_T6bbWW_x09_EM_MAPS/MA5_EM_T6bbWW_x09_GtGrid_SR_8ij80_2ibjet.dat','txt')


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("GtGrid_SR_8ej50_2ibjet")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "GtGrid_SR_8ej50_2ibjet", observedN = 44, expectedBG = 50 , bgError = 10)
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="None"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane([[x,y]]*2)
T1tttt.addSource('efficiencyMap',"./orig/GtGrid_SR_8ej50_2ibjet.txt", "txt")
T1tttt.addSource('obsExclusion', "orig/exclusion_T1tttt.txt", "txt")
T1tttt.dataUrl="http://hepdata.cedar.ac.uk/resource/6095/bjetstream/GtGrid_SR_8ej50_2ibjet.txt"
T1tttt.figure = 'Fig (aux). 11a'
T1tttt.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-04/figaux_11a.png'
#+++++++ next txName block ++++++++++++++
T5WW = dataset.addTxName('T5WW')
T5WW.constraint ="[[['jet','jet'],['W']],[['jet','jet'],['W']]]"
T5WW.conditionDescription =None
T5WW.condition =None
T5WW.massConstraint = None
T5WW.source = "SModelS"
#+++++++ next txName block ++++++++++++++
T5WWoff = dataset.addTxName('T5WWoff')
T5WWoff.constraint ="2.23*[[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
T5WWoff.conditionDescription =None
T5WWoff.condition =None
T5WWoff.massConstraint = [['dm >= 0.0','dm <= 76.']]*2
T5WWoff.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T5WW_1 = T5WW.addMassPlane([[x, 0.5*(x+y), y]]*2)
T5WW_1.addSource('efficiencyMap', 'orig/T5WW_x05/MA5_EM_T5WW_x05_8j50,>=2bjets.dat', 'txt')
T5WW_1.addSource('obsExclusion', 'orig/Exclusion_T5WW_x05.txt', 'txt')
T5WW_1.addSource('obsExclusionM1', 'orig/T5WW_Plus.txt', 'txt')
T5WW_1.addSource('obsExclusionP1', 'orig/T5WW_Minus.txt', 'txt')
T5WW_1.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5WW_2 = T5WW.addMassPlane([[x, 0.05*x + 0.95*y, y]]*2)
T5WW_2.addSource('efficiencyMap','orig/T5WW_x005/MA5_EM_T5WW_Glu005Neu095_8j50,>=2bjets.dat','txt')
T5WW_2.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5WW_3 = T5WW.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5WW_3.addSource('efficiencyMap','orig/T5WW_x095/MA5_EM_T5WW_Glu095Neu005_8j50,>=2bjets.dat','txt')
T5WW_3.dataUrl = None
T5WWoff.addMassPlane(T5WW_1)
T5WWoff.addMassPlane(T5WW_2)
T5WWoff.addMassPlane(T5WW_3)
#+++++++ next txName block ++++++++++++++
T5ZZ = dataset.addTxName('T5ZZ')
T5ZZ.constraint ="[[['jet','jet'],['Z']],[['jet','jet'],['Z']]]"
T5ZZ.conditionDescription =None
T5ZZ.condition =None
T5ZZ.massConstraint = None
T5ZZ.source = "SModelS"
#+++++++ next txName block ++++++++++++++
#+++++++ next mass plane block ++++++++++++++
T5ZZ_1 = T5ZZ.addMassPlane([[x,0.5*(x+y),y]]*2)
T5ZZ_1.addSource('efficiencyMap','orig/T5ZZ_x05/MA5_EM_T5ZZ_x05_8j50,>=2bjets.dat', 'txt')
T5ZZ_1.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5ZZ_2 = T5ZZ.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
T5ZZ_2.addSource('efficiencyMap','orig/T5ZZ_x005/MA5_EM_T5ZZ_Glu005Neu095_8j50,>=2bjets.dat', 'txt')
T5ZZ_2.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5ZZ_3 = T5ZZ.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5ZZ_3.addSource('efficiencyMap','orig/T5ZZ_x095/MA5_EM_T5ZZ_Glu095Neu005_8j50,>=2bjets.dat', 'txt')
T5ZZ_3.dataUrl = None
"""
T5ZZoff = dataset.addTxName('T5ZZoff')
T5ZZoff.constraint ="2.1*[[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
T5ZZoff.conditionDescription =None
T5ZZoff.condition =None
T5ZZoff.massConstraint = [['dm >= 0.0','dm <= 86.']]*2
T5ZZoff.source = "SModelS"
T5ZZoff.addMassPlane(T5ZZ_1)
T5ZZoff.addMassPlane(T5ZZ_2)
T5ZZoff.addMassPlane(T5ZZ_3)
"""
#+++++++ next mass plane block ++++++++++++++
fullpath = "./orig/atlas_susy_2013_04_T2tt_1_EM_MAPS/MA5_EM_T2tt_1_8j50,>=2bjets.dat"
T2tt = dataset.addTxName('T2tt')
T2tt.checked ="None"
T2tt.constraint ="[[['t']],[['t']]]"
T2tt.massConstraint = None
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = "SModelS"
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.addSource('efficiencyMap', fullpath , "txt")
T2tt_1.dataUrl=None
#+++++++ next mass plane block ++++++++++++++
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.checked ="None"
T2ttoff.constraint ="[[['W','b']],[['W','b']]]"
T2ttoff.conditionDescription ="None"
T2ttoff.condition ="None"
T2ttoff.source = "SModelS"
T2ttoff_1 = T2ttoff.addMassPlane([[x,y]]*2)
T2ttoff.massConstraint = [['dm <= 169.0']]*2
T2ttoff_1.addSource('efficiencyMap', fullpath , "txt")
T2ttoff_1.dataUrl=None
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="None"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb.addSource('efficiencyMap',"./orig/atlas_susy_2013_04_T1bbbb_1_EM_MAPS/MA5_EM_T1bbbb_1_8j50,>=2bjets.dat", "txt")
T1bbbb.dataUrl= None
T1btbt = dataset.addTxName('T1btbt')
T1btbt.checked ="None"
T1btbt.constraint ="[[['b','t']],[['b','t']]]"
T1btbt.conditionDescription ="None"
T1btbt.condition ="None"
T1btbt.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T1btbt = T1btbt.addMassPlane([[x,y]]*2)
T1btbt.addSource('efficiencyMap',"./orig/atlas_susy_2013_04_T1btbt_1_EM_MAPS/MA5_EM_T1btbt_1_8j50,>=2bjets.dat", "txt")
T1btbt.dataUrl= None
T5 = dataset.addTxName('T5')
T5.constraint ="[[['jet'],['jet']],[['jet'],['jet']]]"
T5.conditionDescription =None
T5.condition =None
T5.massConstraint = None
T5.source = "SModelS"
## T5.massConstraint = [['dm >= 0.0','dm > 76.']]*2
T5.dataUrl = None
#+++++++ next txName block ++++++++++++++
#+++++++ next mass plane block ++++++++++++++
T5_1 = T5.addMassPlane([[x,0.5*(x+y),y]]*2)
T5_1.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x05_EM_MAPS/MA5_EM_T5_x05_GtGrid_SR_8ej50_2ibjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T5_2 = T5.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
T5_2.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x005_EM_MAPS/MA5_EM_T5_x005_GtGrid_SR_8ej50_2ibjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T5_3 = T5.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5_3.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x095_EM_MAPS/MA5_EM_T5_x095_GtGrid_SR_8ej50_2ibjet.dat', 'txt')
T6bbWW = dataset.addTxName('T6bbWW')
T6bbWW.constraint = "[[['b'],['W']],[['b'],['W']]]"
T6bbWW.conditionDescription =None
T6bbWW.condition =None
T6bbWW.massConstraint = None
T6bbWW.source = "SModelS"
T6bbWW.dataUrl = None
T6bbWW_1 = T6bbWW.addMassPlane([[x, 0.1*x + 0.9*y, y]]*2)
T6bbWW_1.addSource('efficiencyMap', 'orig/atlas_susy_2013_04_T6bbWW_x01_EM_MAPS/MA5_EM_T6bbWW_x01_GtGrid_SR_8ej50_2ibjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T6bbWW_2 = T6bbWW.addMassPlane([[x, 0.5*x + 0.5*y, y]]*2)
T6bbWW_2.addSource('efficiencyMap','orig/atlas_susy_2013_04_T6bbWW_x05_EM_MAPS/MA5_EM_T6bbWW_x05_GtGrid_SR_8ej50_2ibjet.dat','txt')
#+++++++ next mass plane block ++++++++++++++
T6bbWW_3 = T6bbWW.addMassPlane([[x, 0.9*x + 0.1*y, y]]*2)
T6bbWW_3.addSource('efficiencyMap','orig/atlas_susy_2013_04_T6bbWW_x09_EM_MAPS/MA5_EM_T6bbWW_x09_GtGrid_SR_8ej50_2ibjet.dat','txt')


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("GtGrid_SR_7ej80_0bjet")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "GtGrid_SR_7ej80_0bjet", observedN = 12, expectedBG = 11.0 , bgError = 2.2)
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="None"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane([[x,y]]*2)
T1tttt.addSource('efficiencyMap',"./orig/GtGrid_SR_7ej80_0bjet.txt", "txt")
T1tttt.addSource('obsExclusion', "orig/exclusion_T1tttt.txt", "txt")
T1tttt.dataUrl="http://hepdata.cedar.ac.uk/resource/6095/bjetstream/GtGrid_SR_7ej80_0bjet.txt"
T1tttt.figure = 'Fig (aux). 11a'
T1tttt.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-04/figaux_11a.png'
#+++++++ next txName block ++++++++++++++
T5WW = dataset.addTxName('T5WW')
T5WW.constraint ="[[['jet','jet'],['W']],[['jet','jet'],['W']]]"
T5WW.conditionDescription =None
T5WW.condition =None
T5WW.massConstraint = None
T5WW.source = "SModelS"
#+++++++ next txName block ++++++++++++++
T5WWoff = dataset.addTxName('T5WWoff')
T5WWoff.constraint ="2.23*[[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
T5WWoff.conditionDescription =None
T5WWoff.condition =None
T5WWoff.massConstraint = [['dm >= 0.0','dm <= 76.']]*2
T5WWoff.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T5WW_1 = T5WW.addMassPlane([[x, 0.5*(x+y), y]]*2)
T5WW_1.addSource('efficiencyMap', 'orig/T5WW_x05/MA5_EM_T5WW_x05_7j80,0bjet.dat', 'txt')
T5WW_1.addSource('obsExclusion', 'orig/Exclusion_T5WW_x05.txt', 'txt')
T5WW_1.addSource('obsExclusionM1', 'orig/T5WW_Plus.txt', 'txt')
T5WW_1.addSource('obsExclusionP1', 'orig/T5WW_Minus.txt', 'txt')
T5WW_1.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5WW_2 = T5WW.addMassPlane([[x, 0.05*x + 0.95*y, y]]*2)
T5WW_2.addSource('efficiencyMap','orig/T5WW_x005/MA5_EM_T5WW_Glu005Neu095_7j80,0bjet.dat','txt')
T5WW_2.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5WW_3 = T5WW.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5WW_3.addSource('efficiencyMap','orig/T5WW_x095/MA5_EM_T5WW_Glu095Neu005_7j80,0bjet.dat','txt')
T5WW_3.dataUrl = None
T5WWoff.addMassPlane(T5WW_1)
T5WWoff.addMassPlane(T5WW_2)
T5WWoff.addMassPlane(T5WW_3)
#+++++++ next txName block ++++++++++++++
T5ZZ = dataset.addTxName('T5ZZ')
T5ZZ.constraint ="[[['jet','jet'],['Z']],[['jet','jet'],['Z']]]"
T5ZZ.conditionDescription =None
T5ZZ.condition =None
T5ZZ.massConstraint = None
T5ZZ.source = "SModelS"
#+++++++ next txName block ++++++++++++++
#+++++++ next mass plane block ++++++++++++++
T5ZZ_1 = T5ZZ.addMassPlane([[x,0.5*(x+y),y]]*2)
T5ZZ_1.addSource('efficiencyMap','orig/T5ZZ_x05/MA5_EM_T5ZZ_x05_7j80,0bjet.dat', 'txt')
T5ZZ_1.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5ZZ_2 = T5ZZ.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
T5ZZ_2.addSource('efficiencyMap','orig/T5ZZ_x005/MA5_EM_T5ZZ_Glu005Neu095_7j80,0bjet.dat', 'txt')
T5ZZ_2.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5ZZ_3 = T5ZZ.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5ZZ_3.addSource('efficiencyMap','orig/T5ZZ_x095/MA5_EM_T5ZZ_Glu095Neu005_7j80,0bjet.dat', 'txt')
T5ZZ_3.dataUrl = None
"""
T5ZZoff = dataset.addTxName('T5ZZoff')
T5ZZoff.constraint ="2.1*[[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
T5ZZoff.conditionDescription =None
T5ZZoff.condition =None
T5ZZoff.massConstraint = [['dm >= 0.0','dm <= 86.']]*2
T5ZZoff.source = "SModelS"
T5ZZoff.addMassPlane(T5ZZ_1)
T5ZZoff.addMassPlane(T5ZZ_2)
T5ZZoff.addMassPlane(T5ZZ_3)
"""
#+++++++ next mass plane block ++++++++++++++
fullpath = "./orig/atlas_susy_2013_04_T2tt_1_EM_MAPS/MA5_EM_T2tt_1_7j80,0bjet.dat"
T2tt = dataset.addTxName('T2tt')
T2tt.checked ="None"
T2tt.constraint ="[[['t']],[['t']]]"
T2tt.massConstraint = None
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = "SModelS"
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.addSource('efficiencyMap', fullpath , "txt")
T2tt_1.dataUrl=None
#+++++++ next mass plane block ++++++++++++++
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.checked ="None"
T2ttoff.constraint ="[[['W','b']],[['W','b']]]"
T2ttoff.conditionDescription ="None"
T2ttoff.condition ="None"
T2ttoff.source = "SModelS"
T2ttoff_1 = T2ttoff.addMassPlane([[x,y]]*2)
T2ttoff.massConstraint = [['dm <= 169.0']]*2
T2ttoff_1.addSource('efficiencyMap', fullpath , "txt")
T2ttoff_1.dataUrl=None
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="None"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb.addSource('efficiencyMap',"./orig/atlas_susy_2013_04_T1bbbb_1_EM_MAPS/MA5_EM_T1bbbb_1_7j80,0bjet.dat", "txt")
T1bbbb.dataUrl= None
T1btbt = dataset.addTxName('T1btbt')
T1btbt.checked ="None"
T1btbt.constraint ="[[['b','t']],[['b','t']]]"
T1btbt.conditionDescription ="None"
T1btbt.condition ="None"
T1btbt.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T1btbt = T1btbt.addMassPlane([[x,y]]*2)
T1btbt.addSource('efficiencyMap',"./orig/atlas_susy_2013_04_T1btbt_1_EM_MAPS/MA5_EM_T1btbt_1_7j80,0bjet.dat", "txt")
T1btbt.dataUrl= None
T5 = dataset.addTxName('T5')
T5.constraint ="[[['jet'],['jet']],[['jet'],['jet']]]"
T5.conditionDescription =None
T5.condition =None
T5.massConstraint = None
T5.source = "SModelS"
## T5.massConstraint = [['dm >= 0.0','dm > 76.']]*2
T5.dataUrl = None
#+++++++ next txName block ++++++++++++++
#+++++++ next mass plane block ++++++++++++++
T5_1 = T5.addMassPlane([[x,0.5*(x+y),y]]*2)
T5_1.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x05_EM_MAPS/MA5_EM_T5_x05_GtGrid_SR_7ej80_0bjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T5_2 = T5.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
T5_2.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x005_EM_MAPS/MA5_EM_T5_x005_GtGrid_SR_7ej80_0bjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T5_3 = T5.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5_3.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x095_EM_MAPS/MA5_EM_T5_x095_GtGrid_SR_7ej80_0bjet.dat', 'txt')
T6bbWW = dataset.addTxName('T6bbWW')
T6bbWW.constraint = "[[['b'],['W']],[['b'],['W']]]"
T6bbWW.conditionDescription =None
T6bbWW.condition =None
T6bbWW.massConstraint = None
T6bbWW.source = "SModelS"
T6bbWW.dataUrl = None
T6bbWW_1 = T6bbWW.addMassPlane([[x, 0.1*x + 0.9*y, y]]*2)
T6bbWW_1.addSource('efficiencyMap', 'orig/atlas_susy_2013_04_T6bbWW_x01_EM_MAPS/MA5_EM_T6bbWW_x01_GtGrid_SR_7ej80_0bjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T6bbWW_2 = T6bbWW.addMassPlane([[x, 0.5*x + 0.5*y, y]]*2)
T6bbWW_2.addSource('efficiencyMap','orig/atlas_susy_2013_04_T6bbWW_x05_EM_MAPS/MA5_EM_T6bbWW_x05_GtGrid_SR_7ej80_0bjet.dat','txt')
#+++++++ next mass plane block ++++++++++++++
T6bbWW_3 = T6bbWW.addMassPlane([[x, 0.9*x + 0.1*y, y]]*2)
T6bbWW_3.addSource('efficiencyMap','orig/atlas_susy_2013_04_T6bbWW_x09_EM_MAPS/MA5_EM_T6bbWW_x09_GtGrid_SR_7ej80_0bjet.dat','txt')


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("GtGrid_SR_7ej80_1bjet")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "GtGrid_SR_7ej80_1bjet", observedN = 17, expectedBG = 17 , bgError = 6)
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="None"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane([[x,y]]*2)
T1tttt.addSource('efficiencyMap',"./orig/GtGrid_SR_7ej80_1bjet.txt", "txt")
T1tttt.addSource('obsExclusion', "orig/exclusion_T1tttt.txt", "txt")
T1tttt.dataUrl="http://hepdata.cedar.ac.uk/resource/6095/bjetstream/GtGrid_SR_7ej80_1bjet.txt"
T1tttt.figure = 'Fig (aux). 11a'
T1tttt.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-04/figaux_11a.png'
#+++++++ next txName block ++++++++++++++
T5WW = dataset.addTxName('T5WW')
T5WW.constraint ="[[['jet','jet'],['W']],[['jet','jet'],['W']]]"
T5WW.conditionDescription =None
T5WW.condition =None
T5WW.massConstraint = None
T5WW.source = "SModelS"
#+++++++ next txName block ++++++++++++++
T5WWoff = dataset.addTxName('T5WWoff')
T5WWoff.constraint ="2.23*[[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
T5WWoff.conditionDescription =None
T5WWoff.condition =None
T5WWoff.massConstraint = [['dm >= 0.0','dm <= 76.']]*2
T5WWoff.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T5WW_1 = T5WW.addMassPlane([[x, 0.5*(x+y), y]]*2)
T5WW_1.addSource('efficiencyMap', 'orig/T5WW_x05/MA5_EM_T5WW_x05_7j80,1bjet.dat', 'txt')
T5WW_1.addSource('obsExclusion', 'orig/Exclusion_T5WW_x05.txt', 'txt')
T5WW_1.addSource('obsExclusionM1', 'orig/T5WW_Plus.txt', 'txt')
T5WW_1.addSource('obsExclusionP1', 'orig/T5WW_Minus.txt', 'txt')
T5WW_1.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5WW_2 = T5WW.addMassPlane([[x, 0.05*x + 0.95*y, y]]*2)
T5WW_2.addSource('efficiencyMap','orig/T5WW_x005/MA5_EM_T5WW_Glu005Neu095_7j80,1bjet.dat','txt')
T5WW_2.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5WW_3 = T5WW.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5WW_3.addSource('efficiencyMap','orig/T5WW_x095/MA5_EM_T5WW_Glu095Neu005_7j80,1bjet.dat','txt')
T5WW_3.dataUrl = None
T5WWoff.addMassPlane(T5WW_1)
T5WWoff.addMassPlane(T5WW_2)
T5WWoff.addMassPlane(T5WW_3)
#+++++++ next txName block ++++++++++++++
T5ZZ = dataset.addTxName('T5ZZ')
T5ZZ.constraint ="[[['jet','jet'],['Z']],[['jet','jet'],['Z']]]"
T5ZZ.conditionDescription =None
T5ZZ.condition =None
T5ZZ.massConstraint = None
T5ZZ.source = "SModelS"
#+++++++ next txName block ++++++++++++++
#+++++++ next mass plane block ++++++++++++++
T5ZZ_1 = T5ZZ.addMassPlane([[x,0.5*(x+y),y]]*2)
T5ZZ_1.addSource('efficiencyMap','orig/T5ZZ_x05/MA5_EM_T5ZZ_x05_7j80,1bjet.dat', 'txt')
T5ZZ_1.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5ZZ_2 = T5ZZ.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
T5ZZ_2.addSource('efficiencyMap','orig/T5ZZ_x005/MA5_EM_T5ZZ_Glu005Neu095_7j80,1bjet.dat', 'txt')
T5ZZ_2.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5ZZ_3 = T5ZZ.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5ZZ_3.addSource('efficiencyMap','orig/T5ZZ_x095/MA5_EM_T5ZZ_Glu095Neu005_7j80,1bjet.dat', 'txt')
T5ZZ_3.dataUrl = None
"""
T5ZZoff = dataset.addTxName('T5ZZoff')
T5ZZoff.constraint ="2.1*[[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
T5ZZoff.conditionDescription =None
T5ZZoff.condition =None
T5ZZoff.massConstraint = [['dm >= 0.0','dm <= 86.']]*2
T5ZZoff.source = "SModelS"
T5ZZoff.addMassPlane(T5ZZ_1)
T5ZZoff.addMassPlane(T5ZZ_2)
T5ZZoff.addMassPlane(T5ZZ_3)
"""
#+++++++ next mass plane block ++++++++++++++
fullpath = "./orig/atlas_susy_2013_04_T2tt_1_EM_MAPS/MA5_EM_T2tt_1_7j80,1bjet.dat"
T2tt = dataset.addTxName('T2tt')
T2tt.checked ="None"
T2tt.constraint ="[[['t']],[['t']]]"
T2tt.massConstraint = None
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = "SModelS"
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.addSource('efficiencyMap', fullpath , "txt")
T2tt_1.dataUrl=None
#+++++++ next mass plane block ++++++++++++++
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.checked ="None"
T2ttoff.constraint ="[[['W','b']],[['W','b']]]"
T2ttoff.conditionDescription ="None"
T2ttoff.condition ="None"
T2ttoff.source = "SModelS"
T2ttoff_1 = T2ttoff.addMassPlane([[x,y]]*2)
T2ttoff.massConstraint = [['dm <= 169.0']]*2
T2ttoff_1.addSource('efficiencyMap', fullpath , "txt")
T2ttoff_1.dataUrl=None
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="None"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb.addSource('efficiencyMap',"./orig/atlas_susy_2013_04_T1bbbb_1_EM_MAPS/MA5_EM_T1bbbb_1_7j80,1bjet.dat", "txt")
T1bbbb.dataUrl= None
T1btbt = dataset.addTxName('T1btbt')
T1btbt.checked ="None"
T1btbt.constraint ="[[['b','t']],[['b','t']]]"
T1btbt.conditionDescription ="None"
T1btbt.condition ="None"
T1btbt.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T1btbt = T1btbt.addMassPlane([[x,y]]*2)
T1btbt.addSource('efficiencyMap',"./orig/atlas_susy_2013_04_T1btbt_1_EM_MAPS/MA5_EM_T1btbt_1_7j80,1bjet.dat", "txt")
T1btbt.dataUrl= None
T5 = dataset.addTxName('T5')
T5.constraint ="[[['jet'],['jet']],[['jet'],['jet']]]"
T5.conditionDescription =None
T5.condition =None
T5.massConstraint = None
T5.source = "SModelS"
## T5.massConstraint = [['dm >= 0.0','dm > 76.']]*2
T5.dataUrl = None
#+++++++ next txName block ++++++++++++++
#+++++++ next mass plane block ++++++++++++++
T5_1 = T5.addMassPlane([[x,0.5*(x+y),y]]*2)
T5_1.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x05_EM_MAPS/MA5_EM_T5_x05_GtGrid_SR_7ej80_1bjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T5_2 = T5.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
T5_2.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x005_EM_MAPS/MA5_EM_T5_x005_GtGrid_SR_7ej80_1bjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T5_3 = T5.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5_3.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x095_EM_MAPS/MA5_EM_T5_x095_GtGrid_SR_7ej80_1bjet.dat', 'txt')
T6bbWW = dataset.addTxName('T6bbWW')
T6bbWW.constraint = "[[['b'],['W']],[['b'],['W']]]"
T6bbWW.conditionDescription =None
T6bbWW.condition =None
T6bbWW.massConstraint = None
T6bbWW.source = "SModelS"
T6bbWW.dataUrl = None
T6bbWW_1 = T6bbWW.addMassPlane([[x, 0.1*x + 0.9*y, y]]*2)
T6bbWW_1.addSource('efficiencyMap', 'orig/atlas_susy_2013_04_T6bbWW_x01_EM_MAPS/MA5_EM_T6bbWW_x01_GtGrid_SR_7ej80_1bjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T6bbWW_2 = T6bbWW.addMassPlane([[x, 0.5*x + 0.5*y, y]]*2)
T6bbWW_2.addSource('efficiencyMap','orig/atlas_susy_2013_04_T6bbWW_x05_EM_MAPS/MA5_EM_T6bbWW_x05_GtGrid_SR_7ej80_1bjet.dat','txt')
#+++++++ next mass plane block ++++++++++++++
T6bbWW_3 = T6bbWW.addMassPlane([[x, 0.9*x + 0.1*y, y]]*2)
T6bbWW_3.addSource('efficiencyMap','orig/atlas_susy_2013_04_T6bbWW_x09_EM_MAPS/MA5_EM_T6bbWW_x09_GtGrid_SR_7ej80_1bjet.dat','txt')


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("GtGrid_SR_9ej50_2ibjet")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "GtGrid_SR_9ej50_2ibjet", observedN = 7, expectedBG = 8.0 , bgError = 2.7)
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="None"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane([[x,y]]*2)
T1tttt.addSource('efficiencyMap',"./orig/GtGrid_SR_9ej50_2ibjet.txt", "txt")
T1tttt.addSource('obsExclusion', "orig/exclusion_T1tttt.txt", "txt")
T1tttt.dataUrl="http://hepdata.cedar.ac.uk/resource/6095/bjetstream/GtGrid_SR_9ej50_2ibjet.txt"
T1tttt.figure = 'Fig (aux). 11a'
T1tttt.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-04/figaux_11a.png'
#+++++++ next txName block ++++++++++++++
T5WW = dataset.addTxName('T5WW')
T5WW.constraint ="[[['jet','jet'],['W']],[['jet','jet'],['W']]]"
T5WW.conditionDescription =None
T5WW.condition =None
T5WW.massConstraint = None
T5WW.source = "SModelS"
#+++++++ next txName block ++++++++++++++
T5WWoff = dataset.addTxName('T5WWoff')
T5WWoff.constraint ="2.23*[[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
T5WWoff.conditionDescription =None
T5WWoff.condition =None
T5WWoff.massConstraint = [['dm >= 0.0','dm <= 76.']]*2
T5WWoff.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T5WW_1 = T5WW.addMassPlane([[x, 0.5*(x+y), y]]*2)
T5WW_1.addSource('efficiencyMap', 'orig/T5WW_x05/MA5_EM_T5WW_x05_9j50,>=2bjets.dat', 'txt')
T5WW_1.addSource('obsExclusion', 'orig/Exclusion_T5WW_x05.txt', 'txt')
T5WW_1.addSource('obsExclusionM1', 'orig/T5WW_Plus.txt', 'txt')
T5WW_1.addSource('obsExclusionP1', 'orig/T5WW_Minus.txt', 'txt')
T5WW_1.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5WW_2 = T5WW.addMassPlane([[x, 0.05*x + 0.95*y, y]]*2)
T5WW_2.addSource('efficiencyMap','orig/T5WW_x005/MA5_EM_T5WW_Glu005Neu095_9j50,>=2bjets.dat','txt')
T5WW_2.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5WW_3 = T5WW.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5WW_3.addSource('efficiencyMap','orig/T5WW_x095/MA5_EM_T5WW_Glu095Neu005_9j50,>=2bjets.dat','txt')
T5WW_3.dataUrl = None
T5WWoff.addMassPlane(T5WW_1)
T5WWoff.addMassPlane(T5WW_2)
T5WWoff.addMassPlane(T5WW_3)
#+++++++ next txName block ++++++++++++++
T5ZZ = dataset.addTxName('T5ZZ')
T5ZZ.constraint ="[[['jet','jet'],['Z']],[['jet','jet'],['Z']]]"
T5ZZ.conditionDescription =None
T5ZZ.condition =None
T5ZZ.massConstraint = None
T5ZZ.source = "SModelS"
#+++++++ next txName block ++++++++++++++
#+++++++ next mass plane block ++++++++++++++
T5ZZ_1 = T5ZZ.addMassPlane([[x,0.5*(x+y),y]]*2)
T5ZZ_1.addSource('efficiencyMap','orig/T5ZZ_x05/MA5_EM_T5ZZ_x05_9j50,>=2bjets.dat', 'txt')
T5ZZ_1.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5ZZ_2 = T5ZZ.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
T5ZZ_2.addSource('efficiencyMap','orig/T5ZZ_x005/MA5_EM_T5ZZ_Glu005Neu095_9j50,>=2bjets.dat', 'txt')
T5ZZ_2.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5ZZ_3 = T5ZZ.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5ZZ_3.addSource('efficiencyMap','orig/T5ZZ_x095/MA5_EM_T5ZZ_Glu095Neu005_9j50,>=2bjets.dat', 'txt')
T5ZZ_3.dataUrl = None
"""
T5ZZoff = dataset.addTxName('T5ZZoff')
T5ZZoff.constraint ="2.1*[[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
T5ZZoff.conditionDescription =None
T5ZZoff.condition =None
T5ZZoff.massConstraint = [['dm >= 0.0','dm <= 86.']]*2
T5ZZoff.source = "SModelS"
T5ZZoff.addMassPlane(T5ZZ_1)
T5ZZoff.addMassPlane(T5ZZ_2)
T5ZZoff.addMassPlane(T5ZZ_3)
"""
#+++++++ next mass plane block ++++++++++++++
fullpath = "./orig/atlas_susy_2013_04_T2tt_1_EM_MAPS/MA5_EM_T2tt_1_9j50,>=2bjets.dat"
T2tt = dataset.addTxName('T2tt')
T2tt.checked ="None"
T2tt.constraint ="[[['t']],[['t']]]"
T2tt.massConstraint = None
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = "SModelS"
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.addSource('efficiencyMap', fullpath , "txt")
T2tt_1.dataUrl=None
#+++++++ next mass plane block ++++++++++++++
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.checked ="None"
T2ttoff.constraint ="[[['W','b']],[['W','b']]]"
T2ttoff.conditionDescription ="None"
T2ttoff.condition ="None"
T2ttoff.source = "SModelS"
T2ttoff_1 = T2ttoff.addMassPlane([[x,y]]*2)
T2ttoff.massConstraint = [['dm <= 169.0']]*2
T2ttoff_1.addSource('efficiencyMap', fullpath , "txt")
T2ttoff_1.dataUrl=None
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="None"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb.addSource('efficiencyMap',"./orig/atlas_susy_2013_04_T1bbbb_1_EM_MAPS/MA5_EM_T1bbbb_1_9j50,>=2bjets.dat", "txt")
T1bbbb.dataUrl= None
T1btbt = dataset.addTxName('T1btbt')
T1btbt.checked ="None"
T1btbt.constraint ="[[['b','t']],[['b','t']]]"
T1btbt.conditionDescription ="None"
T1btbt.condition ="None"
T1btbt.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T1btbt = T1btbt.addMassPlane([[x,y]]*2)
T1btbt.addSource('efficiencyMap',"./orig/atlas_susy_2013_04_T1btbt_1_EM_MAPS/MA5_EM_T1btbt_1_9j50,>=2bjets.dat", "txt")
T1btbt.dataUrl= None
T5 = dataset.addTxName('T5')
T5.constraint ="[[['jet'],['jet']],[['jet'],['jet']]]"
T5.conditionDescription =None
T5.condition =None
T5.massConstraint = None
T5.source = "SModelS"
## T5.massConstraint = [['dm >= 0.0','dm > 76.']]*2
T5.dataUrl = None
#+++++++ next txName block ++++++++++++++
#+++++++ next mass plane block ++++++++++++++
T5_1 = T5.addMassPlane([[x,0.5*(x+y),y]]*2)
T5_1.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x05_EM_MAPS/MA5_EM_T5_x05_GtGrid_SR_9ej50_2ibjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T5_2 = T5.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
T5_2.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x005_EM_MAPS/MA5_EM_T5_x005_GtGrid_SR_9ej50_2ibjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T5_3 = T5.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5_3.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x095_EM_MAPS/MA5_EM_T5_x095_GtGrid_SR_9ej50_2ibjet.dat', 'txt')
T6bbWW = dataset.addTxName('T6bbWW')
T6bbWW.constraint = "[[['b'],['W']],[['b'],['W']]]"
T6bbWW.conditionDescription =None
T6bbWW.condition =None
T6bbWW.massConstraint = None
T6bbWW.source = "SModelS"
T6bbWW.dataUrl = None
T6bbWW_1 = T6bbWW.addMassPlane([[x, 0.1*x + 0.9*y, y]]*2)
T6bbWW_1.addSource('efficiencyMap', 'orig/atlas_susy_2013_04_T6bbWW_x01_EM_MAPS/MA5_EM_T6bbWW_x01_GtGrid_SR_9ej50_2ibjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T6bbWW_2 = T6bbWW.addMassPlane([[x, 0.5*x + 0.5*y, y]]*2)
T6bbWW_2.addSource('efficiencyMap','orig/atlas_susy_2013_04_T6bbWW_x05_EM_MAPS/MA5_EM_T6bbWW_x05_GtGrid_SR_9ej50_2ibjet.dat','txt')
#+++++++ next mass plane block ++++++++++++++
T6bbWW_3 = T6bbWW.addMassPlane([[x, 0.9*x + 0.1*y, y]]*2)
T6bbWW_3.addSource('efficiencyMap','orig/atlas_susy_2013_04_T6bbWW_x09_EM_MAPS/MA5_EM_T6bbWW_x09_GtGrid_SR_9ej50_2ibjet.dat','txt')


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("GtGrid_SR_7ej80_2ibjet")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "GtGrid_SR_7ej80_2ibjet", observedN = 13, expectedBG = 25 , bgError = 10)
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="None"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane([[x,y]]*2)
T1tttt.addSource('efficiencyMap',"./orig/GtGrid_SR_7ej80_2ibjet.txt", "txt")
T1tttt.addSource('obsExclusion', "orig/exclusion_T1tttt.txt", "txt")
T1tttt.dataUrl="http://hepdata.cedar.ac.uk/resource/6095/bjetstream/GtGrid_SR_7ej80_2ibjet.txt"
T1tttt.figure = 'Fig (aux). 11a'
T1tttt.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-04/figaux_11a.png'
#+++++++ next txName block ++++++++++++++
T5WW = dataset.addTxName('T5WW')
T5WW.constraint ="[[['jet','jet'],['W']],[['jet','jet'],['W']]]"
T5WW.conditionDescription =None
T5WW.condition =None
T5WW.massConstraint = None
T5WW.source = "SModelS"
#+++++++ next txName block ++++++++++++++
T5WWoff = dataset.addTxName('T5WWoff')
T5WWoff.constraint ="2.23*[[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
T5WWoff.conditionDescription =None
T5WWoff.condition =None
T5WWoff.massConstraint = [['dm >= 0.0','dm <= 76.']]*2
T5WWoff.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T5WW_1 = T5WW.addMassPlane([[x, 0.5*(x+y), y]]*2)
T5WW_1.addSource('efficiencyMap', 'orig/T5WW_x05/MA5_EM_T5WW_x05_7j80,>=2bjets.dat', 'txt')
T5WW_1.addSource('obsExclusion', 'orig/Exclusion_T5WW_x05.txt', 'txt')
T5WW_1.addSource('obsExclusionM1', 'orig/T5WW_Plus.txt', 'txt')
T5WW_1.addSource('obsExclusionP1', 'orig/T5WW_Minus.txt', 'txt')
T5WW_1.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5WW_2 = T5WW.addMassPlane([[x, 0.05*x + 0.95*y, y]]*2)
T5WW_2.addSource('efficiencyMap','orig/T5WW_x005/MA5_EM_T5WW_Glu005Neu095_7j80,>=2bjets.dat','txt')
T5WW_2.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5WW_3 = T5WW.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5WW_3.addSource('efficiencyMap','orig/T5WW_x095/MA5_EM_T5WW_Glu095Neu005_7j80,>=2bjets.dat','txt')
T5WW_3.dataUrl = None
T5WWoff.addMassPlane(T5WW_1)
T5WWoff.addMassPlane(T5WW_2)
T5WWoff.addMassPlane(T5WW_3)
#+++++++ next txName block ++++++++++++++
T5ZZ = dataset.addTxName('T5ZZ')
T5ZZ.constraint ="[[['jet','jet'],['Z']],[['jet','jet'],['Z']]]"
T5ZZ.conditionDescription =None
T5ZZ.condition =None
T5ZZ.massConstraint = None
T5ZZ.source = "SModelS"
#+++++++ next txName block ++++++++++++++
#+++++++ next mass plane block ++++++++++++++
T5ZZ_1 = T5ZZ.addMassPlane([[x,0.5*(x+y),y]]*2)
T5ZZ_1.addSource('efficiencyMap','orig/T5ZZ_x05/MA5_EM_T5ZZ_x05_7j80,>=2bjets.dat', 'txt')
T5ZZ_1.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5ZZ_2 = T5ZZ.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
T5ZZ_2.addSource('efficiencyMap','orig/T5ZZ_x005/MA5_EM_T5ZZ_Glu005Neu095_7j80,>=2bjets.dat', 'txt')
T5ZZ_2.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5ZZ_3 = T5ZZ.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5ZZ_3.addSource('efficiencyMap','orig/T5ZZ_x095/MA5_EM_T5ZZ_Glu095Neu005_7j80,>=2bjets.dat', 'txt')
T5ZZ_3.dataUrl = None
"""
T5ZZoff = dataset.addTxName('T5ZZoff')
T5ZZoff.constraint ="2.1*[[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
T5ZZoff.conditionDescription =None
T5ZZoff.condition =None
T5ZZoff.massConstraint = [['dm >= 0.0','dm <= 86.']]*2
T5ZZoff.source = "SModelS"
T5ZZoff.addMassPlane(T5ZZ_1)
T5ZZoff.addMassPlane(T5ZZ_2)
T5ZZoff.addMassPlane(T5ZZ_3)
"""
#+++++++ next mass plane block ++++++++++++++
fullpath = "./orig/atlas_susy_2013_04_T2tt_1_EM_MAPS/MA5_EM_T2tt_1_7j80,>=2bjets.dat"
T2tt = dataset.addTxName('T2tt')
T2tt.checked ="None"
T2tt.constraint ="[[['t']],[['t']]]"
T2tt.massConstraint = None
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = "SModelS"
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.addSource('efficiencyMap', fullpath , "txt")
T2tt_1.dataUrl=None
#+++++++ next mass plane block ++++++++++++++
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.checked ="None"
T2ttoff.constraint ="[[['W','b']],[['W','b']]]"
T2ttoff.conditionDescription ="None"
T2ttoff.condition ="None"
T2ttoff.source = "SModelS"
T2ttoff_1 = T2ttoff.addMassPlane([[x,y]]*2)
T2ttoff.massConstraint = [['dm <= 169.0']]*2
T2ttoff_1.addSource('efficiencyMap', fullpath , "txt")
T2ttoff_1.dataUrl=None
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="None"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb.addSource('efficiencyMap',"./orig/atlas_susy_2013_04_T1bbbb_1_EM_MAPS/MA5_EM_T1bbbb_1_7j80,>=2bjets.dat", "txt")
T1bbbb.dataUrl= None
T1btbt = dataset.addTxName('T1btbt')
T1btbt.checked ="None"
T1btbt.constraint ="[[['b','t']],[['b','t']]]"
T1btbt.conditionDescription ="None"
T1btbt.condition ="None"
T1btbt.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T1btbt = T1btbt.addMassPlane([[x,y]]*2)
T1btbt.addSource('efficiencyMap',"./orig/atlas_susy_2013_04_T1btbt_1_EM_MAPS/MA5_EM_T1btbt_1_7j80,>=2bjets.dat", "txt")
T1btbt.dataUrl= None
T5 = dataset.addTxName('T5')
T5.constraint ="[[['jet'],['jet']],[['jet'],['jet']]]"
T5.conditionDescription =None
T5.condition =None
T5.massConstraint = None
T5.source = "SModelS"
## T5.massConstraint = [['dm >= 0.0','dm > 76.']]*2
T5.dataUrl = None
#+++++++ next txName block ++++++++++++++
#+++++++ next mass plane block ++++++++++++++
T5_1 = T5.addMassPlane([[x,0.5*(x+y),y]]*2)
T5_1.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x05_EM_MAPS/MA5_EM_T5_x05_GtGrid_SR_7ej80_2ibjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T5_2 = T5.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
T5_2.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x005_EM_MAPS/MA5_EM_T5_x005_GtGrid_SR_7ej80_2ibjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T5_3 = T5.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5_3.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x095_EM_MAPS/MA5_EM_T5_x095_GtGrid_SR_7ej80_2ibjet.dat', 'txt')
T6bbWW = dataset.addTxName('T6bbWW')
T6bbWW.constraint = "[[['b'],['W']],[['b'],['W']]]"
T6bbWW.conditionDescription =None
T6bbWW.condition =None
T6bbWW.massConstraint = None
T6bbWW.source = "SModelS"
T6bbWW.dataUrl = None
T6bbWW_1 = T6bbWW.addMassPlane([[x, 0.1*x + 0.9*y, y]]*2)
T6bbWW_1.addSource('efficiencyMap', 'orig/atlas_susy_2013_04_T6bbWW_x01_EM_MAPS/MA5_EM_T6bbWW_x01_GtGrid_SR_7ej80_2ibjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T6bbWW_2 = T6bbWW.addMassPlane([[x, 0.5*x + 0.5*y, y]]*2)
T6bbWW_2.addSource('efficiencyMap','orig/atlas_susy_2013_04_T6bbWW_x05_EM_MAPS/MA5_EM_T6bbWW_x05_GtGrid_SR_7ej80_2ibjet.dat','txt')
#+++++++ next mass plane block ++++++++++++++
T6bbWW_3 = T6bbWW.addMassPlane([[x, 0.9*x + 0.1*y, y]]*2)
T6bbWW_3.addSource('efficiencyMap','orig/atlas_susy_2013_04_T6bbWW_x09_EM_MAPS/MA5_EM_T6bbWW_x09_GtGrid_SR_7ej80_2ibjet.dat','txt')


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("GtGrid_SR_8ej50_1bjet")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "GtGrid_SR_8ej50_1bjet", observedN = 44, expectedBG = 40 , bgError = 10)
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="None"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane([[x,y]]*2)
T1tttt.addSource('efficiencyMap',"./orig/GtGrid_SR_8ej50_1bjet.txt", "txt")
T1tttt.addSource('obsExclusion', "orig/exclusion_T1tttt.txt", "txt")
T1tttt.dataUrl="http://hepdata.cedar.ac.uk/resource/6095/bjetstream/GtGrid_SR_8ej50_1bjet.txt"
T1tttt.figure = 'Fig (aux). 11a'
T1tttt.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-04/figaux_11a.png'
#+++++++ next txName block ++++++++++++++
T5WW = dataset.addTxName('T5WW')
T5WW.constraint ="[[['jet','jet'],['W']],[['jet','jet'],['W']]]"
T5WW.conditionDescription =None
T5WW.condition =None
T5WW.massConstraint = None
T5WW.source = "SModelS"
#+++++++ next txName block ++++++++++++++
T5WWoff = dataset.addTxName('T5WWoff')
T5WWoff.constraint ="2.23*[[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
T5WWoff.conditionDescription =None
T5WWoff.condition =None
T5WWoff.massConstraint = [['dm >= 0.0','dm <= 76.']]*2
T5WWoff.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T5WW_1 = T5WW.addMassPlane([[x, 0.5*(x+y), y]]*2)
T5WW_1.addSource('efficiencyMap', 'orig/T5WW_x05/MA5_EM_T5WW_x05_8j50,1bjet.dat', 'txt')
T5WW_1.addSource('obsExclusion', 'orig/Exclusion_T5WW_x05.txt', 'txt')
T5WW_1.addSource('obsExclusionM1', 'orig/T5WW_Plus.txt', 'txt')
T5WW_1.addSource('obsExclusionP1', 'orig/T5WW_Minus.txt', 'txt')
T5WW_1.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5WW_2 = T5WW.addMassPlane([[x, 0.05*x + 0.95*y, y]]*2)
T5WW_2.addSource('efficiencyMap','orig/T5WW_x005/MA5_EM_T5WW_Glu005Neu095_8j50,1bjet.dat','txt')
T5WW_2.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5WW_3 = T5WW.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5WW_3.addSource('efficiencyMap','orig/T5WW_x095/MA5_EM_T5WW_Glu095Neu005_8j50,1bjet.dat','txt')
T5WW_3.dataUrl = None
T5WWoff.addMassPlane(T5WW_1)
T5WWoff.addMassPlane(T5WW_2)
T5WWoff.addMassPlane(T5WW_3)
#+++++++ next txName block ++++++++++++++
T5ZZ = dataset.addTxName('T5ZZ')
T5ZZ.constraint ="[[['jet','jet'],['Z']],[['jet','jet'],['Z']]]"
T5ZZ.conditionDescription =None
T5ZZ.condition =None
T5ZZ.massConstraint = None
T5ZZ.source = "SModelS"
#+++++++ next txName block ++++++++++++++
#+++++++ next mass plane block ++++++++++++++
T5ZZ_1 = T5ZZ.addMassPlane([[x,0.5*(x+y),y]]*2)
T5ZZ_1.addSource('efficiencyMap','orig/T5ZZ_x05/MA5_EM_T5ZZ_x05_8j50,1bjet.dat', 'txt')
T5ZZ_1.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5ZZ_2 = T5ZZ.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
T5ZZ_2.addSource('efficiencyMap','orig/T5ZZ_x005/MA5_EM_T5ZZ_Glu005Neu095_8j50,1bjet.dat', 'txt')
T5ZZ_2.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5ZZ_3 = T5ZZ.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5ZZ_3.addSource('efficiencyMap','orig/T5ZZ_x095/MA5_EM_T5ZZ_Glu095Neu005_8j50,1bjet.dat', 'txt')
T5ZZ_3.dataUrl = None
"""
T5ZZoff = dataset.addTxName('T5ZZoff')
T5ZZoff.constraint ="2.1*[[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
T5ZZoff.conditionDescription =None
T5ZZoff.condition =None
T5ZZoff.massConstraint = [['dm >= 0.0','dm <= 86.']]*2
T5ZZoff.source = "SModelS"
T5ZZoff.addMassPlane(T5ZZ_1)
T5ZZoff.addMassPlane(T5ZZ_2)
T5ZZoff.addMassPlane(T5ZZ_3)
"""
#+++++++ next mass plane block ++++++++++++++
fullpath = "./orig/atlas_susy_2013_04_T2tt_1_EM_MAPS/MA5_EM_T2tt_1_8j50,1bjet.dat"
T2tt = dataset.addTxName('T2tt')
T2tt.checked ="None"
T2tt.constraint ="[[['t']],[['t']]]"
T2tt.massConstraint = None
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = "SModelS"
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.addSource('efficiencyMap', fullpath , "txt")
T2tt_1.dataUrl=None
#+++++++ next mass plane block ++++++++++++++
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.checked ="None"
T2ttoff.constraint ="[[['W','b']],[['W','b']]]"
T2ttoff.conditionDescription ="None"
T2ttoff.condition ="None"
T2ttoff.source = "SModelS"
T2ttoff_1 = T2ttoff.addMassPlane([[x,y]]*2)
T2ttoff.massConstraint = [['dm <= 169.0']]*2
T2ttoff_1.addSource('efficiencyMap', fullpath , "txt")
T2ttoff_1.dataUrl=None
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="None"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb.addSource('efficiencyMap',"./orig/atlas_susy_2013_04_T1bbbb_1_EM_MAPS/MA5_EM_T1bbbb_1_8j50,1bjet.dat", "txt")
T1bbbb.dataUrl= None
T1btbt = dataset.addTxName('T1btbt')
T1btbt.checked ="None"
T1btbt.constraint ="[[['b','t']],[['b','t']]]"
T1btbt.conditionDescription ="None"
T1btbt.condition ="None"
T1btbt.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T1btbt = T1btbt.addMassPlane([[x,y]]*2)
T1btbt.addSource('efficiencyMap',"./orig/atlas_susy_2013_04_T1btbt_1_EM_MAPS/MA5_EM_T1btbt_1_8j50,1bjet.dat", "txt")
T1btbt.dataUrl= None
T5 = dataset.addTxName('T5')
T5.constraint ="[[['jet'],['jet']],[['jet'],['jet']]]"
T5.conditionDescription =None
T5.condition =None
T5.massConstraint = None
T5.source = "SModelS"
## T5.massConstraint = [['dm >= 0.0','dm > 76.']]*2
T5.dataUrl = None
#+++++++ next txName block ++++++++++++++
#+++++++ next mass plane block ++++++++++++++
T5_1 = T5.addMassPlane([[x,0.5*(x+y),y]]*2)
T5_1.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x05_EM_MAPS/MA5_EM_T5_x05_GtGrid_SR_8ej50_1bjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T5_2 = T5.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
T5_2.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x005_EM_MAPS/MA5_EM_T5_x005_GtGrid_SR_8ej50_1bjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T5_3 = T5.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5_3.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x095_EM_MAPS/MA5_EM_T5_x095_GtGrid_SR_8ej50_1bjet.dat', 'txt')
T6bbWW = dataset.addTxName('T6bbWW')
T6bbWW.constraint = "[[['b'],['W']],[['b'],['W']]]"
T6bbWW.conditionDescription =None
T6bbWW.condition =None
T6bbWW.massConstraint = None
T6bbWW.source = "SModelS"
T6bbWW.dataUrl = None
T6bbWW_1 = T6bbWW.addMassPlane([[x, 0.1*x + 0.9*y, y]]*2)
T6bbWW_1.addSource('efficiencyMap', 'orig/atlas_susy_2013_04_T6bbWW_x01_EM_MAPS/MA5_EM_T6bbWW_x01_GtGrid_SR_8ej50_1bjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T6bbWW_2 = T6bbWW.addMassPlane([[x, 0.5*x + 0.5*y, y]]*2)
T6bbWW_2.addSource('efficiencyMap','orig/atlas_susy_2013_04_T6bbWW_x05_EM_MAPS/MA5_EM_T6bbWW_x05_GtGrid_SR_8ej50_1bjet.dat','txt')
#+++++++ next mass plane block ++++++++++++++
T6bbWW_3 = T6bbWW.addMassPlane([[x, 0.9*x + 0.1*y, y]]*2)
T6bbWW_3.addSource('efficiencyMap','orig/atlas_susy_2013_04_T6bbWW_x09_EM_MAPS/MA5_EM_T6bbWW_x09_GtGrid_SR_8ej50_1bjet.dat','txt')


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("GtGrid_SR_8ij80_1bjet")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "GtGrid_SR_8ij80_1bjet", observedN = 1, expectedBG = 1.5 , bgError = 0.9)
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="None"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane([[x,y]]*2)
T1tttt.addSource('efficiencyMap',"./orig/GtGrid_SR_8ij80_1bjet.txt", "txt")
T1tttt.addSource('obsExclusion', "orig/exclusion_T1tttt.txt", "txt")
T1tttt.dataUrl="http://hepdata.cedar.ac.uk/resource/6095/bjetstream/GtGrid_SR_8ij80_1bjet.txt"
T1tttt.figure = 'Fig (aux). 11a'
T1tttt.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-04/figaux_11a.png'
#+++++++ next txName block ++++++++++++++
T5WW = dataset.addTxName('T5WW')
T5WW.constraint ="[[['jet','jet'],['W']],[['jet','jet'],['W']]]"
T5WW.conditionDescription =None
T5WW.condition =None
T5WW.massConstraint = None
T5WW.source = "SModelS"
#+++++++ next txName block ++++++++++++++
T5WWoff = dataset.addTxName('T5WWoff')
T5WWoff.constraint ="2.23*[[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
T5WWoff.conditionDescription =None
T5WWoff.condition =None
T5WWoff.massConstraint = [['dm >= 0.0','dm <= 76.']]*2
T5WWoff.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T5WW_1 = T5WW.addMassPlane([[x, 0.5*(x+y), y]]*2)
T5WW_1.addSource('efficiencyMap', 'orig/T5WW_x05/MA5_EM_T5WW_x05_>=8j80,1bjet.dat', 'txt')
T5WW_1.addSource('obsExclusion', 'orig/Exclusion_T5WW_x05.txt', 'txt')
T5WW_1.addSource('obsExclusionM1', 'orig/T5WW_Plus.txt', 'txt')
T5WW_1.addSource('obsExclusionP1', 'orig/T5WW_Minus.txt', 'txt')
T5WW_1.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5WW_2 = T5WW.addMassPlane([[x, 0.05*x + 0.95*y, y]]*2)
T5WW_2.addSource('efficiencyMap','orig/T5WW_x005/MA5_EM_T5WW_Glu005Neu095_>=8j80,1bjet.dat','txt')
T5WW_2.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5WW_3 = T5WW.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5WW_3.addSource('efficiencyMap','orig/T5WW_x095/MA5_EM_T5WW_Glu095Neu005_>=8j80,1bjet.dat','txt')
T5WW_3.dataUrl = None
T5WWoff.addMassPlane(T5WW_1)
T5WWoff.addMassPlane(T5WW_2)
T5WWoff.addMassPlane(T5WW_3)
#+++++++ next txName block ++++++++++++++
T5ZZ = dataset.addTxName('T5ZZ')
T5ZZ.constraint ="[[['jet','jet'],['Z']],[['jet','jet'],['Z']]]"
T5ZZ.conditionDescription =None
T5ZZ.condition =None
T5ZZ.massConstraint = None
T5ZZ.source = "SModelS"
#+++++++ next txName block ++++++++++++++
#+++++++ next mass plane block ++++++++++++++
T5ZZ_1 = T5ZZ.addMassPlane([[x,0.5*(x+y),y]]*2)
T5ZZ_1.addSource('efficiencyMap','orig/T5ZZ_x05/MA5_EM_T5ZZ_x05_>=8j80,1bjet.dat', 'txt')
T5ZZ_1.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5ZZ_2 = T5ZZ.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
T5ZZ_2.addSource('efficiencyMap','orig/T5ZZ_x005/MA5_EM_T5ZZ_Glu005Neu095_>=8j80,1bjet.dat', 'txt')
T5ZZ_2.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5ZZ_3 = T5ZZ.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5ZZ_3.addSource('efficiencyMap','orig/T5ZZ_x095/MA5_EM_T5ZZ_Glu095Neu005_>=8j80,1bjet.dat', 'txt')
T5ZZ_3.dataUrl = None
"""
T5ZZoff = dataset.addTxName('T5ZZoff')
T5ZZoff.constraint ="2.1*[[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
T5ZZoff.conditionDescription =None
T5ZZoff.condition =None
T5ZZoff.massConstraint = [['dm >= 0.0','dm <= 86.']]*2
T5ZZoff.source = "SModelS"
T5ZZoff.addMassPlane(T5ZZ_1)
T5ZZoff.addMassPlane(T5ZZ_2)
T5ZZoff.addMassPlane(T5ZZ_3)
"""
#+++++++ next mass plane block ++++++++++++++
fullpath = "./orig/atlas_susy_2013_04_T2tt_1_EM_MAPS/MA5_EM_T2tt_1_>=8j80,1bjet.dat"
T2tt = dataset.addTxName('T2tt')
T2tt.checked ="None"
T2tt.constraint ="[[['t']],[['t']]]"
T2tt.massConstraint = None
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = "SModelS"
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.addSource('efficiencyMap', fullpath , "txt")
T2tt_1.dataUrl=None
#+++++++ next mass plane block ++++++++++++++
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.checked ="None"
T2ttoff.constraint ="[[['W','b']],[['W','b']]]"
T2ttoff.conditionDescription ="None"
T2ttoff.condition ="None"
T2ttoff.source = "SModelS"
T2ttoff_1 = T2ttoff.addMassPlane([[x,y]]*2)
T2ttoff.massConstraint = [['dm <= 169.0']]*2
T2ttoff_1.addSource('efficiencyMap', fullpath , "txt")
T2ttoff_1.dataUrl=None
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="None"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb.addSource('efficiencyMap',"./orig/atlas_susy_2013_04_T1bbbb_1_EM_MAPS/MA5_EM_T1bbbb_1_>=8j80,1bjet.dat", "txt")
T1bbbb.dataUrl= None
T1btbt = dataset.addTxName('T1btbt')
T1btbt.checked ="None"
T1btbt.constraint ="[[['b','t']],[['b','t']]]"
T1btbt.conditionDescription ="None"
T1btbt.condition ="None"
T1btbt.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T1btbt = T1btbt.addMassPlane([[x,y]]*2)
T1btbt.addSource('efficiencyMap',"./orig/atlas_susy_2013_04_T1btbt_1_EM_MAPS/MA5_EM_T1btbt_1_>=8j80,1bjet.dat", "txt")
T1btbt.dataUrl= None
T5 = dataset.addTxName('T5')
T5.constraint ="[[['jet'],['jet']],[['jet'],['jet']]]"
T5.conditionDescription =None
T5.condition =None
T5.massConstraint = None
T5.source = "SModelS"
## T5.massConstraint = [['dm >= 0.0','dm > 76.']]*2
T5.dataUrl = None
#+++++++ next txName block ++++++++++++++
#+++++++ next mass plane block ++++++++++++++
T5_1 = T5.addMassPlane([[x,0.5*(x+y),y]]*2)
T5_1.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x05_EM_MAPS/MA5_EM_T5_x05_GtGrid_SR_8ij80_1bjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T5_2 = T5.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
T5_2.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x005_EM_MAPS/MA5_EM_T5_x005_GtGrid_SR_8ij80_1bjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T5_3 = T5.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5_3.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x095_EM_MAPS/MA5_EM_T5_x095_GtGrid_SR_8ij80_1bjet.dat', 'txt')
T6bbWW = dataset.addTxName('T6bbWW')
T6bbWW.constraint = "[[['b'],['W']],[['b'],['W']]]"
T6bbWW.conditionDescription =None
T6bbWW.condition =None
T6bbWW.massConstraint = None
T6bbWW.source = "SModelS"
T6bbWW.dataUrl = None
T6bbWW_1 = T6bbWW.addMassPlane([[x, 0.1*x + 0.9*y, y]]*2)
T6bbWW_1.addSource('efficiencyMap', 'orig/atlas_susy_2013_04_T6bbWW_x01_EM_MAPS/MA5_EM_T6bbWW_x01_GtGrid_SR_8ij80_1bjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T6bbWW_2 = T6bbWW.addMassPlane([[x, 0.5*x + 0.5*y, y]]*2)
T6bbWW_2.addSource('efficiencyMap','orig/atlas_susy_2013_04_T6bbWW_x05_EM_MAPS/MA5_EM_T6bbWW_x05_GtGrid_SR_8ij80_1bjet.dat','txt')
#+++++++ next mass plane block ++++++++++++++
T6bbWW_3 = T6bbWW.addMassPlane([[x, 0.9*x + 0.1*y, y]]*2)
T6bbWW_3.addSource('efficiencyMap','orig/atlas_susy_2013_04_T6bbWW_x09_EM_MAPS/MA5_EM_T6bbWW_x09_GtGrid_SR_8ij80_1bjet.dat','txt')


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("GtGrid_SR_9ej50_1bjet")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "GtGrid_SR_9ej50_1bjet", observedN = 8, expectedBG = 6.1 , bgError = 1.7)
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="None"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane([[x,y]]*2)
T1tttt.addSource('efficiencyMap',"./orig/GtGrid_SR_9ej50_1bjet.txt", "txt")
T1tttt.addSource('obsExclusion', "orig/exclusion_T1tttt.txt", "txt")
T1tttt.dataUrl="http://hepdata.cedar.ac.uk/resource/6095/bjetstream/GtGrid_SR_9ej50_1bjet.txt"
T1tttt.figure = 'Fig (aux). 11a'
T1tttt.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-04/figaux_11a.png'
#+++++++ next txName block ++++++++++++++
T5WW = dataset.addTxName('T5WW')
T5WW.constraint ="[[['jet','jet'],['W']],[['jet','jet'],['W']]]"
T5WW.conditionDescription =None
T5WW.condition =None
T5WW.massConstraint = None
T5WW.source = "SModelS"
#+++++++ next txName block ++++++++++++++
T5WWoff = dataset.addTxName('T5WWoff')
T5WWoff.constraint ="2.23*[[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
T5WWoff.conditionDescription =None
T5WWoff.condition =None
T5WWoff.massConstraint = [['dm >= 0.0','dm <= 76.']]*2
T5WWoff.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T5WW_1 = T5WW.addMassPlane([[x, 0.5*(x+y), y]]*2)
T5WW_1.addSource('efficiencyMap', 'orig/T5WW_x05/MA5_EM_T5WW_x05_9j50,1bjet.dat', 'txt')
T5WW_1.addSource('obsExclusion', 'orig/Exclusion_T5WW_x05.txt', 'txt')
T5WW_1.addSource('obsExclusionM1', 'orig/T5WW_Plus.txt', 'txt')
T5WW_1.addSource('obsExclusionP1', 'orig/T5WW_Minus.txt', 'txt')
T5WW_1.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5WW_2 = T5WW.addMassPlane([[x, 0.05*x + 0.95*y, y]]*2)
T5WW_2.addSource('efficiencyMap','orig/T5WW_x005/MA5_EM_T5WW_Glu005Neu095_9j50,1bjet.dat','txt')
T5WW_2.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5WW_3 = T5WW.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5WW_3.addSource('efficiencyMap','orig/T5WW_x095/MA5_EM_T5WW_Glu095Neu005_9j50,1bjet.dat','txt')
T5WW_3.dataUrl = None
T5WWoff.addMassPlane(T5WW_1)
T5WWoff.addMassPlane(T5WW_2)
T5WWoff.addMassPlane(T5WW_3)
#+++++++ next txName block ++++++++++++++
T5ZZ = dataset.addTxName('T5ZZ')
T5ZZ.constraint ="[[['jet','jet'],['Z']],[['jet','jet'],['Z']]]"
T5ZZ.conditionDescription =None
T5ZZ.condition =None
T5ZZ.massConstraint = None
T5ZZ.source = "SModelS"
#+++++++ next txName block ++++++++++++++
#+++++++ next mass plane block ++++++++++++++
T5ZZ_1 = T5ZZ.addMassPlane([[x,0.5*(x+y),y]]*2)
T5ZZ_1.addSource('efficiencyMap','orig/T5ZZ_x05/MA5_EM_T5ZZ_x05_9j50,1bjet.dat', 'txt')
T5ZZ_1.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5ZZ_2 = T5ZZ.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
T5ZZ_2.addSource('efficiencyMap','orig/T5ZZ_x005/MA5_EM_T5ZZ_Glu005Neu095_9j50,1bjet.dat', 'txt')
T5ZZ_2.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5ZZ_3 = T5ZZ.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5ZZ_3.addSource('efficiencyMap','orig/T5ZZ_x095/MA5_EM_T5ZZ_Glu095Neu005_9j50,1bjet.dat', 'txt')
T5ZZ_3.dataUrl = None
"""
T5ZZoff = dataset.addTxName('T5ZZoff')
T5ZZoff.constraint ="2.1*[[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
T5ZZoff.conditionDescription =None
T5ZZoff.condition =None
T5ZZoff.massConstraint = [['dm >= 0.0','dm <= 86.']]*2
T5ZZoff.source = "SModelS"
T5ZZoff.addMassPlane(T5ZZ_1)
T5ZZoff.addMassPlane(T5ZZ_2)
T5ZZoff.addMassPlane(T5ZZ_3)
"""
#+++++++ next mass plane block ++++++++++++++
fullpath = "./orig/atlas_susy_2013_04_T2tt_1_EM_MAPS/MA5_EM_T2tt_1_9j50,1bjet.dat"
T2tt = dataset.addTxName('T2tt')
T2tt.checked ="None"
T2tt.constraint ="[[['t']],[['t']]]"
T2tt.massConstraint = None
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = "SModelS"
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.addSource('efficiencyMap', fullpath , "txt")
T2tt_1.dataUrl=None
#+++++++ next mass plane block ++++++++++++++
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.checked ="None"
T2ttoff.constraint ="[[['W','b']],[['W','b']]]"
T2ttoff.conditionDescription ="None"
T2ttoff.condition ="None"
T2ttoff.source = "SModelS"
T2ttoff_1 = T2ttoff.addMassPlane([[x,y]]*2)
T2ttoff.massConstraint = [['dm <= 169.0']]*2
T2ttoff_1.addSource('efficiencyMap', fullpath , "txt")
T2ttoff_1.dataUrl=None
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="None"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb.addSource('efficiencyMap',"./orig/atlas_susy_2013_04_T1bbbb_1_EM_MAPS/MA5_EM_T1bbbb_1_9j50,1bjet.dat", "txt")
T1bbbb.dataUrl= None
T1btbt = dataset.addTxName('T1btbt')
T1btbt.checked ="None"
T1btbt.constraint ="[[['b','t']],[['b','t']]]"
T1btbt.conditionDescription ="None"
T1btbt.condition ="None"
T1btbt.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T1btbt = T1btbt.addMassPlane([[x,y]]*2)
T1btbt.addSource('efficiencyMap',"./orig/atlas_susy_2013_04_T1btbt_1_EM_MAPS/MA5_EM_T1btbt_1_9j50,1bjet.dat", "txt")
T1btbt.dataUrl= None
T5 = dataset.addTxName('T5')
T5.constraint ="[[['jet'],['jet']],[['jet'],['jet']]]"
T5.conditionDescription =None
T5.condition =None
T5.massConstraint = None
T5.source = "SModelS"
## T5.massConstraint = [['dm >= 0.0','dm > 76.']]*2
T5.dataUrl = None
#+++++++ next txName block ++++++++++++++
#+++++++ next mass plane block ++++++++++++++
T5_1 = T5.addMassPlane([[x,0.5*(x+y),y]]*2)
T5_1.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x05_EM_MAPS/MA5_EM_T5_x05_GtGrid_SR_9ej50_1bjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T5_2 = T5.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
T5_2.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x005_EM_MAPS/MA5_EM_T5_x005_GtGrid_SR_9ej50_1bjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T5_3 = T5.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5_3.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x095_EM_MAPS/MA5_EM_T5_x095_GtGrid_SR_9ej50_1bjet.dat', 'txt')
T6bbWW = dataset.addTxName('T6bbWW')
T6bbWW.constraint = "[[['b'],['W']],[['b'],['W']]]"
T6bbWW.conditionDescription =None
T6bbWW.condition =None
T6bbWW.massConstraint = None
T6bbWW.source = "SModelS"
T6bbWW.dataUrl = None
T6bbWW_1 = T6bbWW.addMassPlane([[x, 0.1*x + 0.9*y, y]]*2)
T6bbWW_1.addSource('efficiencyMap', 'orig/atlas_susy_2013_04_T6bbWW_x01_EM_MAPS/MA5_EM_T6bbWW_x01_GtGrid_SR_9ej50_1bjet.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T6bbWW_2 = T6bbWW.addMassPlane([[x, 0.5*x + 0.5*y, y]]*2)
T6bbWW_2.addSource('efficiencyMap','orig/atlas_susy_2013_04_T6bbWW_x05_EM_MAPS/MA5_EM_T6bbWW_x05_GtGrid_SR_9ej50_1bjet.dat','txt')
#+++++++ next mass plane block ++++++++++++++
T6bbWW_3 = T6bbWW.addMassPlane([[x, 0.9*x + 0.1*y, y]]*2)
T6bbWW_3.addSource('efficiencyMap','orig/atlas_susy_2013_04_T6bbWW_x09_EM_MAPS/MA5_EM_T6bbWW_x09_GtGrid_SR_9ej50_1bjet.dat','txt')


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("GtGrid_SR_10ij50_bjetblind")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "GtGrid_SR_10ij50_bjetblind", observedN = 3, expectedBG = 1.37 , bgError = 0.35)
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="None"
T1tttt.constraint ="[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane([[x,y]]*2)
T1tttt.addSource('efficiencyMap',"./orig/GtGrid_SR_10ij50_bjetblind.txt", "txt")
T1tttt.addSource('obsExclusion', "orig/exclusion_T1tttt.txt", "txt")
T1tttt.dataUrl="http://hepdata.cedar.ac.uk/resource/6095/ins1247060_resources.tar.gz"
T1tttt.figure = 'Fig (aux). 11a'
T1tttt.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-04/figaux_11a.png'
#+++++++ next txName block ++++++++++++++
T5WW = dataset.addTxName('T5WW')
T5WW.constraint ="[[['jet','jet'],['W']],[['jet','jet'],['W']]]"
T5WW.conditionDescription =None
T5WW.condition =None
T5WW.massConstraint = None
T5WW.source = "SModelS"
#+++++++ next txName block ++++++++++++++
T5WWoff = dataset.addTxName('T5WWoff')
T5WWoff.constraint ="2.23*[[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
T5WWoff.conditionDescription =None
T5WWoff.condition =None
T5WWoff.massConstraint = [['dm >= 0.0','dm <= 76.']]*2
T5WWoff.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T5WW_1 = T5WW.addMassPlane([[x, 0.5*(x+y), y]]*2)
T5WW_1.addSource('efficiencyMap', 'orig/T5WW_x05/MA5_EM_T5WW_x05_>=10j50.dat', 'txt')
T5WW_1.addSource('obsExclusion', 'orig/Exclusion_T5WW_x05.txt', 'txt')
T5WW_1.addSource('obsExclusionM1', 'orig/T5WW_Plus.txt', 'txt')
T5WW_1.addSource('obsExclusionP1', 'orig/T5WW_Minus.txt', 'txt')
T5WW_1.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5WW_2 = T5WW.addMassPlane([[x, 0.05*x + 0.95*y, y]]*2)
T5WW_2.addSource('efficiencyMap','orig/T5WW_x005/MA5_EM_T5WW_Glu005Neu095_>=10j50.dat','txt')
T5WW_2.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5WW_3 = T5WW.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5WW_3.addSource('efficiencyMap','orig/T5WW_x095/MA5_EM_T5WW_Glu095Neu005_>=10j50.dat','txt')
T5WW_3.dataUrl = None
T5WWoff.addMassPlane(T5WW_1)
T5WWoff.addMassPlane(T5WW_2)
T5WWoff.addMassPlane(T5WW_3)
#+++++++ next txName block ++++++++++++++
T5ZZ = dataset.addTxName('T5ZZ')
T5ZZ.constraint ="[[['jet','jet'],['Z']],[['jet','jet'],['Z']]]"
T5ZZ.conditionDescription =None
T5ZZ.condition =None
T5ZZ.massConstraint = None
T5ZZ.source = "SModelS"
#+++++++ next txName block ++++++++++++++
#+++++++ next mass plane block ++++++++++++++
T5ZZ_1 = T5ZZ.addMassPlane([[x,0.5*(x+y),y]]*2)
T5ZZ_1.addSource('efficiencyMap','orig/T5ZZ_x05/MA5_EM_T5ZZ_x05_>=10j50.dat', 'txt')
T5ZZ_1.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5ZZ_2 = T5ZZ.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
T5ZZ_2.addSource('efficiencyMap','orig/T5ZZ_x005/MA5_EM_T5ZZ_Glu005Neu095_>=10j50.dat', 'txt')
T5ZZ_2.dataUrl = None
#+++++++ next mass plane block ++++++++++++++
T5ZZ_3 = T5ZZ.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5ZZ_3.addSource('efficiencyMap','orig/T5ZZ_x095/MA5_EM_T5ZZ_Glu095Neu005_>=10j50.dat', 'txt')
T5ZZ_3.dataUrl = None
"""
T5ZZoff = dataset.addTxName('T5ZZoff')
T5ZZoff.constraint ="2.1*[[['jet','jet'],['jet','jet']],[['jet','jet'],['jet','jet']]]"
T5ZZoff.conditionDescription =None
T5ZZoff.condition =None
T5ZZoff.massConstraint = [['dm >= 0.0','dm <= 86.']]*2
T5ZZoff.source = "SModelS"
T5ZZoff.addMassPlane(T5ZZ_1)
T5ZZoff.addMassPlane(T5ZZ_2)
T5ZZoff.addMassPlane(T5ZZ_3)
"""
#+++++++ next mass plane block ++++++++++++++
fullpath = "./orig/atlas_susy_2013_04_T2tt_1_EM_MAPS/MA5_EM_T2tt_1_>=10j50.dat"
T2tt = dataset.addTxName('T2tt')
T2tt.checked ="None"
T2tt.constraint ="[[['t']],[['t']]]"
T2tt.massConstraint = None
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = "SModelS"
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.addSource('efficiencyMap', fullpath , "txt")
T2tt_1.dataUrl=None
#+++++++ next mass plane block ++++++++++++++
T2ttoff = dataset.addTxName('T2ttoff')
T2ttoff.checked ="None"
T2ttoff.constraint ="[[['W','b']],[['W','b']]]"
T2ttoff.conditionDescription ="None"
T2ttoff.condition ="None"
T2ttoff.source = "SModelS"
T2ttoff_1 = T2ttoff.addMassPlane([[x,y]]*2)
T2ttoff.massConstraint = [['dm <= 169.0']]*2
T2ttoff_1.addSource('efficiencyMap', fullpath , "txt")
T2ttoff_1.dataUrl=None
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="None"
T1bbbb.constraint ="[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane([[x,y]]*2)
T1bbbb.addSource('efficiencyMap',"./orig/atlas_susy_2013_04_T1bbbb_1_EM_MAPS/MA5_EM_T1bbbb_1_>=10j50.dat", "txt")
T1bbbb.dataUrl= None
T1btbt = dataset.addTxName('T1btbt')
T1btbt.checked ="None"
T1btbt.constraint ="[[['b','t']],[['b','t']]]"
T1btbt.conditionDescription ="None"
T1btbt.condition ="None"
T1btbt.source = "SModelS"
#+++++++ next mass plane block ++++++++++++++
T1btbt = T1btbt.addMassPlane([[x,y]]*2)
T1btbt.addSource('efficiencyMap',"./orig/atlas_susy_2013_04_T1btbt_1_EM_MAPS/MA5_EM_T1btbt_1_>=10j50.dat", "txt")
T1btbt.dataUrl= None
T5 = dataset.addTxName('T5')
T5.constraint ="[[['jet'],['jet']],[['jet'],['jet']]]"
T5.conditionDescription =None
T5.condition =None
T5.massConstraint = None
T5.source = "SModelS"
## T5.massConstraint = [['dm >= 0.0','dm > 76.']]*2
T5.dataUrl = None
#+++++++ next txName block ++++++++++++++
#+++++++ next mass plane block ++++++++++++++
T5_1 = T5.addMassPlane([[x,0.5*(x+y),y]]*2)
T5_1.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x05_EM_MAPS/MA5_EM_T5_x05_GtGrid_SR_10ij50_bjetblind.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T5_2 = T5.addMassPlane([[x,0.05*x + 0.95*y,y]]*2)
T5_2.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x005_EM_MAPS/MA5_EM_T5_x005_GtGrid_SR_10ij50_bjetblind.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T5_3 = T5.addMassPlane([[x, 0.95*x + 0.05*y, y]]*2)
T5_3.addSource('efficiencyMap','orig/atlas_susy_2013_04_T5_x095_EM_MAPS/MA5_EM_T5_x095_GtGrid_SR_10ij50_bjetblind.dat', 'txt')
T6bbWW = dataset.addTxName('T6bbWW')
T6bbWW.constraint = "[[['b'],['W']],[['b'],['W']]]"
T6bbWW.conditionDescription =None
T6bbWW.condition =None
T6bbWW.massConstraint = None
T6bbWW.source = "SModelS"
T6bbWW.dataUrl = None
T6bbWW_1 = T6bbWW.addMassPlane([[x, 0.1*x + 0.9*y, y]]*2)
T6bbWW_1.addSource('efficiencyMap', 'orig/atlas_susy_2013_04_T6bbWW_x01_EM_MAPS/MA5_EM_T6bbWW_x01_GtGrid_SR_10ij50_bjetblind.dat', 'txt')
#+++++++ next mass plane block ++++++++++++++
T6bbWW_2 = T6bbWW.addMassPlane([[x, 0.5*x + 0.5*y, y]]*2)
T6bbWW_2.addSource('efficiencyMap','orig/atlas_susy_2013_04_T6bbWW_x05_EM_MAPS/MA5_EM_T6bbWW_x05_GtGrid_SR_10ij50_bjetblind.dat','txt')
#+++++++ next mass plane block ++++++++++++++
T6bbWW_3 = T6bbWW.addMassPlane([[x, 0.9*x + 0.1*y, y]]*2)
T6bbWW_3.addSource('efficiencyMap','orig/atlas_susy_2013_04_T6bbWW_x09_EM_MAPS/MA5_EM_T6bbWW_x09_GtGrid_SR_10ij50_bjetblind.dat','txt')

databaseCreator.create()
