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
info = MetaInfoInput('ATLAS-SUSY-2013-09')
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-09/'
info.sqrts = 8
info.lumi = '20.3'
info.prettyName = '2 SS leptons'
info.private = False
info.arxiv = 'http://arxiv.org/abs/1404.2500'
info.contact = "ATLAS collaboration"
info.publication = 'http://link.springer.com/article/10.1007/JHEP06(2014)035'
info.supersedes = 'ATLAS-CONF-2013-007; ATLAS-CONF-2012-151'



#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR3b")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR3b", observedN = 1, expectedBG = 2.2 , bgError = 0.8, upperLimit = '1.927E-01*fb', expectedUpperLimit = '2.438E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = "None"
T1tttt.condition = "None"
T1tttt.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
# T1tttt_1.figure ="No Figure"
T1tttt_1.figureUrl = None
T1tttt_1.addSource('obsExclusion', "orig/exclusion_T1tttt.txt", "txt")
T1tttt_1.addSource('efficiencyMap', "orig/EffMap_T1tttt_SR3b.txt", "txt")
T1tttt_1.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1289225'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR3Lhigh")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR3Lhigh", observedN = 2, expectedBG = 2.5 , bgError = 0.9, upperLimit = '2.388E-01*fb', expectedUpperLimit = '2.395E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = "None"
T1tttt.condition = "None"
T1tttt.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure ="No Figure"
T1tttt_1.figureUrl ="No Figure"
T1tttt_1.addSource('obsExclusion', "orig/exclusion_T1tttt.txt", "txt")
T1tttt_1.addSource('efficiencyMap', "orig/EffMap_T1tttt_SR3Lhigh.txt", "txt")
T1tttt_1.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1289225'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR3Llow")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR3Llow", observedN = 6, expectedBG = 4.3 , bgError = 2.1, upperLimit = '4.312E-01*fb', expectedUpperLimit = '3.307E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = "None"
T1tttt.condition = "None"
T1tttt.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure ="No Figure"
T1tttt_1.figureUrl ="No Figure"
T1tttt_1.addSource('obsExclusion', "orig/exclusion_T1tttt.txt", "txt")
T1tttt_1.addSource('efficiencyMap', "orig/EffMap_T1tttt_SR3Llow.txt", "txt")
T1tttt_1.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1289225'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR0b")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR0b", observedN = 14, expectedBG = 6.5 , bgError = 2.3, upperLimit = '8.011E-01*fb', expectedUpperLimit = '3.794E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = "None"
T1tttt.condition = "None"
T1tttt.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure ="No Figure"
T1tttt_1.figureUrl ="No Figure"
T1tttt_1.addSource('obsExclusion', "orig/exclusion_T1tttt.txt", "txt")
T1tttt_1.addSource('efficiencyMap', "orig/EffMap_T1tttt_SR0b.txt", "txt")
T1tttt_1.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1289225'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR1b")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR1b", observedN = 10, expectedBG = 4.7 , bgError = 2.1, upperLimit = '6.434E-01*fb', expectedUpperLimit = '3.220E-01*fb')
#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription = "None"
T1tttt.condition = "None"
T1tttt.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1tttt_1 = T1tttt.addMassPlane([[x,y]]*2)
T1tttt_1.figure ="No Figure"
T1tttt_1.figureUrl ="No Figure"
T1tttt_1.addSource('obsExclusion', "orig/exclusion_T1tttt.txt", "txt")
T1tttt_1.addSource('efficiencyMap', "orig/EffMap_T1tttt_SR1b.txt", "txt")
T1tttt_1.dataUrl = 'http://hepdata.cedar.ac.uk/view/ins1289225'

databaseCreator.create()
