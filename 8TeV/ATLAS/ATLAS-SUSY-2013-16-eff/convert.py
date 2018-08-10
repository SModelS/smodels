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
info = MetaInfoInput('ATLAS-SUSY-2013-16')
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-16/'
info.sqrts = '8'
info.lumi = '20.1'
info.prettyName = '0 lepton + 6 (2 b-)jets + Etmiss'
info.private = False
info.arxiv = 'http://arxiv.org/abs/1406.1122'
info.contact = ''
info.publication = 'http://link.springer.com/article/10.1007%2FJHEP09%282014%29015'
info.comment ='EM given for All SRs for one topo, only for 3 for the other. Only these 3 are implemented.'
info.supersedes ='ATLAS-CONF-2013-024'



#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SRA4")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SRA4", observedN = 4, expectedBG = 2.4 , bgError = 0.7, upperLimit = '3.5106E-01*fb', expectedUpperLimit = '2.3921E-01*fb')
#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.constraint = "[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.addSource('obsExclusion', "./orig/T2tt_Exclusion.txt", "txt")
T2tt_1.addSource('efficiencyMap','orig/EffMap_T2tt_SRA4.txt', 'txt')
T2tt_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1299143/d27/input http://hepdata.cedar.ac.uk/view/ins1299143/d27/input"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SRC2")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SRC2", observedN = 30, expectedBG = 34 , bgError = 5, upperLimit = '6.7114E-01*fb', expectedUpperLimit = '8.1205E-01*fb')
#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.constraint = "[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.addSource('obsExclusion', "./orig/T2tt_Exclusion.txt", "txt")
T2tt_1.addSource('efficiencyMap','orig/EffMap_T2tt_SRC2.txt', 'txt')
T2tt_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1299143/d27/input http://hepdata.cedar.ac.uk/view/ins1299143/d27/input"
#+++++++ next txName block ++++++++++++++
T6bbWW = dataset.addTxName('T6bbWW')
T6bbWW.constraint ="[[['b'],['W']],[['b'],['W']]]"
T6bbWW.conditionDescription ="None"
T6bbWW.condition ="None"
T6bbWW.massConstraint = None
T6bbWW.source = 'ATLAS'
#+++++++ next txName block ++++++++++++++
T6bbWWoff = dataset.addTxName('T6bbWWoff')
T6bbWWoff.condition ="None"
T6bbWWoff.constraint ="[[['b'],['jet','jet']],[['b'],['jet','jet']]]"
T6bbWWoff.conditionDescription ="None"
T6bbWWoff.massConstraint = [['dm >= 0.0','dm <= 76.']]*2
T6bbWWoff.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T6bbWW_1 = T6bbWW.addMassPlane([[x,2*y,y]]*2)
T6bbWW_1.addSource('obsExclusion', "./orig/T6bbWW_Exclusion.txt" ,"txt")
T6bbWW_1.addSource('efficiencyMap', "./orig/EffMap_T6bbWW_SRC2.txt", "txt")
T6bbWW_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1299143/next http://hepdata.cedar.ac.uk/view/ins1299143/d29/input"
T6bbWWoff.addMassPlane(T6bbWW_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SRC3")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SRC3", observedN = 15, expectedBG = 20.3 , bgError = 3.0, upperLimit = '4.2061E-01*fb', expectedUpperLimit = '5.8823E-01*fb')
#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.constraint = "[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.addSource('obsExclusion', "./orig/T2tt_Exclusion.txt", "txt")
T2tt_1.addSource('efficiencyMap','orig/EffMap_T2tt_SRC3.txt', 'txt')
T2tt_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1299143/d27/input http://hepdata.cedar.ac.uk/view/ins1299143/d27/input"
#+++++++ next txName block ++++++++++++++
T6bbWW = dataset.addTxName('T6bbWW')
T6bbWW.constraint ="[[['b'],['W']],[['b'],['W']]]"
T6bbWW.conditionDescription ="None"
T6bbWW.condition ="None"
T6bbWW.massConstraint = None
T6bbWW.source = 'ATLAS'
#+++++++ next txName block ++++++++++++++
T6bbWWoff = dataset.addTxName('T6bbWWoff')
T6bbWWoff.condition ="None"
T6bbWWoff.constraint ="[[['b'],['jet','jet']],[['b'],['jet','jet']]]"
T6bbWWoff.conditionDescription ="None"
T6bbWWoff.massConstraint = [['dm >= 0.0','dm <= 76.']]*2
T6bbWWoff.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T6bbWW_1 = T6bbWW.addMassPlane([[x,2*y,y]]*2)
T6bbWW_1.addSource('obsExclusion', "./orig/T6bbWW_Exclusion.txt" ,"txt")
T6bbWW_1.addSource('efficiencyMap', "./orig/EffMap_T6bbWW_SRC3.txt", "txt")
T6bbWW_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1299143/next http://hepdata.cedar.ac.uk/view/ins1299143/d29/input"
T6bbWWoff.addMassPlane(T6bbWW_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SRC1")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SRC1", observedN = 59, expectedBG = 68 , bgError = 7, upperLimit = '8.3525E-01*fb', expectedUpperLimit = '1.1191E+00*fb')
#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.constraint = "[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.addSource('obsExclusion', "./orig/T2tt_Exclusion.txt", "txt")
T2tt_1.addSource('efficiencyMap','orig/EffMap_T2tt_SRC1.txt', 'txt')
T2tt_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1299143/d27/input http://hepdata.cedar.ac.uk/view/ins1299143/d27/input"
#+++++++ next txName block ++++++++++++++
T6bbWW = dataset.addTxName('T6bbWW')
T6bbWW.constraint ="[[['b'],['W']],[['b'],['W']]]"
T6bbWW.conditionDescription ="None"
T6bbWW.condition ="None"
T6bbWW.massConstraint = None
T6bbWW.source = 'ATLAS'
#+++++++ next txName block ++++++++++++++
T6bbWWoff = dataset.addTxName('T6bbWWoff')
T6bbWWoff.condition ="None"
T6bbWWoff.constraint ="[[['b'],['jet','jet']],[['b'],['jet','jet']]]"
T6bbWWoff.conditionDescription ="None"
T6bbWWoff.massConstraint = [['dm >= 0.0','dm <= 76.']]*2
T6bbWWoff.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T6bbWW_1 = T6bbWW.addMassPlane([[x,2*y,y]]*2)
T6bbWW_1.addSource('obsExclusion', "./orig/T6bbWW_Exclusion.txt" ,"txt")
T6bbWW_1.addSource('efficiencyMap', "./orig/EffMap_T6bbWW_SRC1.txt", "txt")
T6bbWW_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1299143/next http://hepdata.cedar.ac.uk/view/ins1299143/d29/input"
T6bbWWoff.addMassPlane(T6bbWW_1)


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SRA1")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SRA1", observedN = 11, expectedBG = 15.8 , bgError = 1.9, upperLimit = '3.4146E-01*fb', expectedUpperLimit = '4.7789E-01*fb')
#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.constraint = "[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.addSource('obsExclusion', "./orig/T2tt_Exclusion.txt", "txt")
T2tt_1.addSource('efficiencyMap','orig/EffMap_T2tt_SRA1.txt', 'txt')
T2tt_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1299143/d27/input http://hepdata.cedar.ac.uk/view/ins1299143/d27/input"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SRA2")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SRA2", observedN = 4, expectedBG = 4.1 , bgError = 0.8, upperLimit = '3.0024E-01*fb', expectedUpperLimit = '2.9905E-01*fb')
#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.constraint = "[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.addSource('obsExclusion', "./orig/T2tt_Exclusion.txt", "txt")
T2tt_1.addSource('efficiencyMap','orig/EffMap_T2tt_SRA2.txt', 'txt')
T2tt_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1299143/d27/input http://hepdata.cedar.ac.uk/view/ins1299143/d27/input"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SRA3")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SRA3", observedN = 5, expectedBG = 4.1 , bgError = 0.9, upperLimit = '3.5508E-01*fb', expectedUpperLimit = '3.0183E-01*fb')
#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.constraint = "[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.addSource('obsExclusion', "./orig/T2tt_Exclusion.txt", "txt")
T2tt_1.addSource('efficiencyMap','orig/EffMap_T2tt_SRA3.txt', 'txt')
T2tt_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1299143/d27/input http://hepdata.cedar.ac.uk/view/ins1299143/d27/input"

databaseCreator.create()
