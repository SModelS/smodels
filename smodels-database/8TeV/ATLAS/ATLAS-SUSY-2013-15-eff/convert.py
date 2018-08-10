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
info = MetaInfoInput('ATLAS-SUSY-2013-15')
info.url = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-15/"
info.sqrts = 8
info.lumi = 20.3
info.prettyName = '1 lepton + 4 (1 b-)jets + Etmiss'
info.private = False
info.arxiv = "http://arxiv.org/abs/1407.0583"
info.contact = "ATLAS Collaboration"
info.publication = "http://link.springer.com/article/10.1007/JHEP11(2014)118"
info.comment = "Only T2tt 'can' be impl., but 1 SR EM is missing (3 are given out of 4). tNdiag is MET>150,mT>140 ."
info.supersedes = 'ATLAS-CONF-2013-037; CONF-2012-166'



#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("tNboost")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "tNboost", observedN = 5, expectedBG = 3.3 , bgError = 0.7, upperLimit = '3.726E-01*fb', expectedUpperLimit = '2.678E-01*fb')
#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.constraint = "[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.addSource('obsExclusion', "orig/exclusion_T2tt.txt", "txt")
T2tt_1.addSource('efficiencyMap','orig/EffMap_T2tt_tNboost.txt', 'txt')
T2tt_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1304456"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("tNdiag")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "tNdiag", observedN = 217, expectedBG = 236 , bgError = 29, upperLimit = '2.636E+00*fb', expectedUpperLimit = '3.223E+00*fb')
#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.constraint = "[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.addSource('obsExclusion', "orig/exclusion_T2tt.txt", "txt")
T2tt_1.addSource('efficiencyMap','orig/EffMap_T2tt_tNdiag.txt', 'txt')
T2tt_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1304456"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("tNmed")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "tNmed", observedN = 12, expectedBG = 13 , bgError = 2.2, upperLimit = '4.463E-01*fb', expectedUpperLimit = '4.844E-01*fb')
#+++++++ next txName block ++++++++++++++
T2tt = dataset.addTxName('T2tt')
T2tt.constraint = "[[['t']],[['t']]]"
T2tt.conditionDescription ="None"
T2tt.condition ="None"
T2tt.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane([[x,y]]*2)
T2tt_1.addSource('obsExclusion', "orig/exclusion_T2tt.txt", "txt")
T2tt_1.addSource('efficiencyMap','orig/EffMap_T2tt_tNmed.txt', 'txt')
T2tt_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1304456"

databaseCreator.create()
