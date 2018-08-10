#!/usr/bin/env python

"""
.. module:: convert
   :synopsis: used to create info.txt and the <txname>.txt files.

"""
import sys
import os
import argparse

argparser = argparse.ArgumentParser(description = \
'create info.txt, txname.txt, twiki.txt and sms.py')
argparser.add_argument ( '-utilsPath', '--utilsPath',
			                   help = 'path to the package smodels_utils',\
					               type = str )
argparser.add_argument ( '-smodelsPath', '--smodelsPath',
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

#+++++++ global info block ++++++++++++++
info = MetaInfoInput('ATLAS-SUSY-2015-06')
info.url = 'http://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2015-06/'
info.sqrts = 13
info.lumi = 3.2
info.prettyName = '2-6 jets, 0 lep'
info.private = False
info.arxiv =  'https://arxiv.org/abs/1605.03814'
info.contact = 'ATLAS collaboration'
info.publication = 'http://link.springer.com/article/10.1140/epjc/s10052-016-4184-8'
info.comment = 'UL analyses does not have digital data available'
#info.supersedes =
#info.supersededBy =



#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR6jm")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR6jm", observedN = 4, expectedBG = 6.9 , bgError = 1.5)
T2 = dataset.addTxName('T2')
T2.checked = 'NO'
T2.constraint ="[[['jet']],[['jet']]]"
T2.conditionDescription = None
T2.condition = None
T2.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T2_1 = T2.addMassPlane( [[x,y]]*2 )
#---- new efficiency map -----
#----figure----
T2_1.figureUrl  = "http://hepdata.cedar.ac.uk/resource/1458270/figAuxiliaryFigure24b.png"
T2_1.figure  = "Fig. Aux. 24b"
#----exclusion source----
T2_1.addSource( 'obsExclusion', 'orig/Obs_Line_T2.dat', 'txt', objectName = None, index = None )
T2_1.addSource( 'expExclusion', 'orig/Exp_Line_T2.dat', 'txt', objectName = None, index = None )
T2_1.addSource( 'efficiencyMap', 'orig/AccXEff_T2_SR6jm.dat', 'txt', objectName = None, index = None )
#T2bb_1.expExclusionM1.setSource( path, type, objectName = None, index = None )
#T2bb_1.expExclusionP1.setSource( path, type, objectName = None, index = None )
T2_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1458270/d75/input"
T1 = dataset.addTxName('T1')
T1.checked = 'NO'
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription = None
T1.condition = None
T1.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T1_1 = T1.addMassPlane( [[x,y]]*2 )
#---- new efficiency map -----
#----figure----
T1_1.figureUrl  = "http://hepdata.cedar.ac.uk/resource/1458270/figAuxiliaryFigure31b.png"
T1_1.figure  = "Fig. Aux. 31b"
#----exclusion source----
T1_1.addSource( 'obsExclusion', 'orig/Obs_Line_T1.dat', 'txt', objectName = None, index = None )
T1_1.addSource( 'expExclusion', 'orig/Exp_Line_T1.dat', 'txt', objectName = None, index = None )
T1_1.addSource( 'efficiencyMap', 'orig/AccXEff_T1_SR6jm.dat', 'txt', objectName = None, index = None )
#T2bb_1.expExclusionM1.setSource( path, type, objectName = None, index = None )
#T2bb_1.expExclusionP1.setSource( path, type, objectName = None, index = None )
T1_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1458270/d89/input"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR2jt")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR2jt", observedN = 26, expectedBG = 23 , bgError = 4)
T2 = dataset.addTxName('T2')
T2.checked = 'NO'
T2.constraint ="[[['jet']],[['jet']]]"
T2.conditionDescription = None
T2.condition = None
T2.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T2_1 = T2.addMassPlane( [[x,y]]*2 )
#---- new efficiency map -----
#----figure----
T2_1.figureUrl  = "http://hepdata.cedar.ac.uk/resource/1458270/figAuxiliaryFigure21b.png"
T2_1.figure  = "Fig. Aux. 21b"
#----exclusion source----
T2_1.addSource( 'obsExclusion', 'orig/Obs_Line_T2.dat', 'txt', objectName = None, index = None )
T2_1.addSource( 'expExclusion', 'orig/Exp_Line_T2.dat', 'txt', objectName = None, index = None )
T2_1.addSource( 'efficiencyMap', 'orig/AccXEff_T2_SR2jt.dat', 'txt', objectName = None, index = None )
#T2bb_1.expExclusionM1.setSource( path, type, objectName = None, index = None )
#T2bb_1.expExclusionP1.setSource( path, type, objectName = None, index = None )
T2_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1458270/d69/input"
T1 = dataset.addTxName('T1')
T1.checked = 'NO'
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription = None
T1.condition = None
T1.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T1_1 = T1.addMassPlane( [[x,y]]*2 )
#---- new efficiency map -----
#----figure----
T1_1.figureUrl  = "http://hepdata.cedar.ac.uk/resource/1458270/figAuxiliaryFigure28b.png"
T1_1.figure  = "Fig. Aux. 28b"
#----exclusion source----
T1_1.addSource( 'obsExclusion', 'orig/Obs_Line_T1.dat', 'txt', objectName = None, index = None )
T1_1.addSource( 'expExclusion', 'orig/Exp_Line_T1.dat', 'txt', objectName = None, index = None )
T1_1.addSource( 'efficiencyMap', 'orig/AccXEff_T1_SR2jt.dat', 'txt', objectName = None, index = None )
#T2bb_1.expExclusionM1.setSource( path, type, objectName = None, index = None )
#T2bb_1.expExclusionP1.setSource( path, type, objectName = None, index = None )
T1_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1458270/d81/input"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR4jt")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR4jt", observedN = 7, expectedBG = 4.1 , bgError = 1.1)
T2 = dataset.addTxName('T2')
T2.checked = 'NO'
T2.constraint ="[[['jet']],[['jet']]]"
T2.conditionDescription = None
T2.condition = None
T2.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T2_1 = T2.addMassPlane( [[x,y]]*2 )
#---- new efficiency map -----
#----figure----
T2_1.figureUrl  = "http://hepdata.cedar.ac.uk/resource/1458270/figAuxiliaryFigure22b.png"
T2_1.figure  = "Fig. Aux. 22b"
#----exclusion source----
T2_1.addSource( 'obsExclusion', 'orig/Obs_Line_T2.dat', 'txt', objectName = None, index = None )
T2_1.addSource( 'expExclusion', 'orig/Exp_Line_T2.dat', 'txt', objectName = None, index = None )
T2_1.addSource( 'efficiencyMap', 'orig/AccXEff_T2_SR4jt.dat', 'txt', objectName = None, index = None )
#T2bb_1.expExclusionM1.setSource( path, type, objectName = None, index = None )
#T2bb_1.expExclusionP1.setSource( path, type, objectName = None, index = None )
T2_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1458270/d71/input"
T1 = dataset.addTxName('T1')
T1.checked = 'NO'
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription = None
T1.condition = None
T1.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T1_1 = T1.addMassPlane( [[x,y]]*2 )
#---- new efficiency map -----
#----figure----
T1_1.figureUrl  = "http://hepdata.cedar.ac.uk/resource/1458270/figAuxiliaryFigure29b.png"
T1_1.figure  = "Fig. Aux. 29b"
#----exclusion source----
T1_1.addSource( 'obsExclusion', 'orig/Obs_Line_T1.dat', 'txt', objectName = None, index = None )
T1_1.addSource( 'expExclusion', 'orig/Exp_Line_T1.dat', 'txt', objectName = None, index = None )
T1_1.addSource( 'efficiencyMap', 'orig/AccXEff_T1_SR4jt.dat', 'txt', objectName = None, index = None )
#T2bb_1.expExclusionM1.setSource( path, type, objectName = None, index = None )
#T2bb_1.expExclusionP1.setSource( path, type, objectName = None, index = None )
T1_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1458270/d85/input"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR2jl")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR2jl", observedN = 263, expectedBG = 283 , bgError = 24)
T2 = dataset.addTxName('T2')
T2.checked = 'NO'
T2.constraint ="[[['jet']],[['jet']]]"
T2.conditionDescription = None
T2.condition = None
T2.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T2_1 = T2.addMassPlane( [[x,y]]*2 )
#---- new efficiency map -----
#----figure----
T2_1.figureUrl  = "http://hepdata.cedar.ac.uk/resource/1458270/figAuxiliaryFigure19b.png"
T2_1.figure  = "Fig. Aux. 19b"
#----exclusion source----
T2_1.addSource( 'obsExclusion', 'orig/Obs_Line_T2.dat', 'txt', objectName = None, index = None )
T2_1.addSource( 'expExclusion', 'orig/Exp_Line_T2.dat', 'txt', objectName = None, index = None )
T2_1.addSource( 'efficiencyMap', 'orig/AccXEff_T2_SR2jl.dat', 'txt', objectName = None, index = None )
#T2bb_1.expExclusionM1.setSource( path, type, objectName = None, index = None )
#T2bb_1.expExclusionP1.setSource( path, type, objectName = None, index = None )
T2_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1458270/d65/input"
T1 = dataset.addTxName('T1')
T1.checked = 'NO'
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription = None
T1.condition = None
T1.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T1_1 = T1.addMassPlane( [[x,y]]*2 )
#---- new efficiency map -----
#----figure----
T1_1.figureUrl  = "http://hepdata.cedar.ac.uk/resource/1458270/figAuxiliaryFigure26b.png"
T1_1.figure  = "Fig. Aux. 26b"
#----exclusion source----
T1_1.addSource( 'obsExclusion', 'orig/Obs_Line_T1.dat', 'txt', objectName = None, index = None )
T1_1.addSource( 'expExclusion', 'orig/Exp_Line_T1.dat', 'txt', objectName = None, index = None )
T1_1.addSource( 'efficiencyMap', 'orig/AccXEff_T1_SR2jl.dat', 'txt', objectName = None, index = None )
#T2bb_1.expExclusionM1.setSource( path, type, objectName = None, index = None )
#T2bb_1.expExclusionP1.setSource( path, type, objectName = None, index = None )
T1_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1458270/d79/input"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR2jm")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR2jm", observedN = 191, expectedBG = 191 , bgError = 21)
T2 = dataset.addTxName('T2')
T2.checked = 'NO'
T2.constraint ="[[['jet']],[['jet']]]"
T2.conditionDescription = None
T2.condition = None
T2.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T2_1 = T2.addMassPlane( [[x,y]]*2 )
#---- new efficiency map -----
#----figure----
T2_1.figureUrl  = "http://hepdata.cedar.ac.uk/resource/1458270/figAuxiliaryFigure20b.png"
T2_1.figure  = "Fig. Aux. 20b"
#----exclusion source----
T2_1.addSource( 'obsExclusion', 'orig/Obs_Line_T2.dat', 'txt', objectName = None, index = None )
T2_1.addSource( 'expExclusion', 'orig/Exp_Line_T2.dat', 'txt', objectName = None, index = None )
T2_1.addSource( 'efficiencyMap', 'orig/AccXEff_T2_SR2jm.dat', 'txt', objectName = None, index = None )
#T2bb_1.expExclusionM1.setSource( path, type, objectName = None, index = None )
#T2bb_1.expExclusionP1.setSource( path, type, objectName = None, index = None )
T2_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1458270/d67/input"
T1 = dataset.addTxName('T1')
T1.checked = 'NO'
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription = None
T1.condition = None
T1.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T1_1 = T1.addMassPlane( [[x,y]]*2 )
#---- new efficiency map -----
#----figure----
T1_1.figureUrl  = "http://hepdata.cedar.ac.uk/resource/1458270/figAuxiliaryFigure27b.png"
T1_1.figure  = "Fig. Aux. 27b"
#----exclusion source----
T1_1.addSource( 'obsExclusion', 'orig/Obs_Line_T1.dat', 'txt', objectName = None, index = None )
T1_1.addSource( 'expExclusion', 'orig/Exp_Line_T1.dat', 'txt', objectName = None, index = None )
T1_1.addSource( 'efficiencyMap', 'orig/AccXEff_T1_SR2jm.dat', 'txt', objectName = None, index = None )
#T2bb_1.expExclusionM1.setSource( path, type, objectName = None, index = None )
#T2bb_1.expExclusionP1.setSource( path, type, objectName = None, index = None )
T1_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1458270/d79/input"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR5j")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR5j", observedN = 7, expectedBG = 13.2 , bgError = 2.2)
T2 = dataset.addTxName('T2')
T2.checked = 'NO'
T2.constraint ="[[['jet']],[['jet']]]"
T2.conditionDescription = None
T2.condition = None
T2.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T2_1 = T2.addMassPlane( [[x,y]]*2 )
#---- new efficiency map -----
#----figure----
T2_1.figureUrl  = "http://hepdata.cedar.ac.uk/resource/1458270/figAuxiliaryFigure23b.png"
T2_1.figure  = "Fig. Aux. 23b"
#----exclusion source----
T2_1.addSource( 'obsExclusion', 'orig/Obs_Line_T2.dat', 'txt', objectName = None, index = None )
T2_1.addSource( 'expExclusion', 'orig/Exp_Line_T2.dat', 'txt', objectName = None, index = None )
T2_1.addSource( 'efficiencyMap', 'orig/AccXEff_T2_SR5j.dat', 'txt', objectName = None, index = None )
#T2bb_1.expExclusionM1.setSource( path, type, objectName = None, index = None )
#T2bb_1.expExclusionP1.setSource( path, type, objectName = None, index = None )
T2_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1458270/d73/input"
T1 = dataset.addTxName('T1')
T1.checked = 'NO'
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription = None
T1.condition = None
T1.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T1_1 = T1.addMassPlane( [[x,y]]*2 )
#---- new efficiency map -----
#----figure----
T1_1.figureUrl  = "http://hepdata.cedar.ac.uk/resource/1458270/figAuxiliaryFigure30b.png"
T1_1.figure  = "Fig. Aux. 30b"
#----exclusion source----
T1_1.addSource( 'obsExclusion', 'orig/Obs_Line_T1.dat', 'txt', objectName = None, index = None )
T1_1.addSource( 'expExclusion', 'orig/Exp_Line_T1.dat', 'txt', objectName = None, index = None )
T1_1.addSource( 'efficiencyMap', 'orig/AccXEff_T1_SR5j.dat', 'txt', objectName = None, index = None )
#T2bb_1.expExclusionM1.setSource( path, type, objectName = None, index = None )
#T2bb_1.expExclusionP1.setSource( path, type, objectName = None, index = None )
T1_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1458270/d87/input"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR6jt")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR6jt", observedN = 3, expectedBG = 4.2 , bgError = 1.2)
T2 = dataset.addTxName('T2')
T2.checked = 'NO'
T2.constraint ="[[['jet']],[['jet']]]"
T2.conditionDescription = None
T2.condition = None
T2.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T2_1 = T2.addMassPlane( [[x,y]]*2 )
#---- new efficiency map -----
#----figure----
T2_1.figureUrl  = "http://hepdata.cedar.ac.uk/resource/1458270/figAuxiliaryFigure25b.png"
T2_1.figure  = "Fig.Aux. 25b"
#----exclusion source----
T2_1.addSource( 'obsExclusion', 'orig/Obs_Line_T2.dat', 'txt', objectName = None, index = None )
T2_1.addSource( 'expExclusion', 'orig/Exp_Line_T2.dat', 'txt', objectName = None, index = None )
T2_1.addSource( 'efficiencyMap', 'orig/AccXEff_T2_SR6jt.dat', 'txt', objectName = None, index = None )
#T2bb_1.expExclusionM1.setSource( path, type, objectName = None, index = None )
#T2bb_1.expExclusionP1.setSource( path, type, objectName = None, index = None )
T2_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1458270/d77/input"
T1 = dataset.addTxName('T1')
T1.checked = 'NO'
T1.constraint ="[[['jet','jet']],[['jet','jet']]]"
T1.conditionDescription = None
T1.condition = None
T1.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T1_1 = T1.addMassPlane( [[x,y]]*2 )
#---- new efficiency map -----
#----figure----
T1_1.figureUrl  = "http://hepdata.cedar.ac.uk/resource/1458270/figAuxiliaryFigure32b.png"
T1_1.figure  = "Fig. Aux. 32b"
#----exclusion source----
T1_1.addSource( 'obsExclusion', 'orig/Obs_Line_T1.dat', 'txt', objectName = None, index = None )
T1_1.addSource( 'expExclusion', 'orig/Exp_Line_T1.dat', 'txt', objectName = None, index = None )
T1_1.addSource( 'efficiencyMap', 'orig/AccXEff_T1_SR6jt.dat', 'txt', objectName = None, index = None )
#T2bb_1.expExclusionM1.setSource( path, type, objectName = None, index = None )
#T2bb_1.expExclusionP1.setSource( path, type, objectName = None, index = None )
T1_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1458270/d91/input"
databaseCreator.create()

