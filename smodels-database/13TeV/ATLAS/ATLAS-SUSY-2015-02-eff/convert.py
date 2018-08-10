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
info = MetaInfoInput('ATLAS-SUSY-2015-02')
info.url = 'http://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2015-02/'
info.sqrts = 13
info.lumi = 3.2
info.prettyName = 'Top 1l'
info.private = False
info.arxiv =  'https://arxiv.org/abs/1606.03903'
info.source = "ATLAS"
info.publication ='http://journals.aps.org/prd/abstract/10.1103/PhysRevD.94.052009' 
info.comment = 'Very weird exclusion due to best CLs selection, and fluctuation in the Obs data' 
#info.supersedes =
#info.supersededBy =



#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR2")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR2", observedN = 1, expectedBG = 1.25 , bgError = 0.26)
T2tt = dataset.addTxName('T2tt')
T2tt.checked = 'NO'
T2tt.constraint ="[[['t']],[['t']]]"
T2tt.conditionDescription = None
T2tt.condition = None
T2tt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane( [[x,y]]*2 )
#---- new efficiency map -----
#----figure----
T2tt_1.figureUrl  = "http://hepdata.cedar.ac.uk/resource/1469069/figAuxiliaryFigure10b.png; http://hepdata.cedar.ac.uk/resource/1469069/figAuxiliaryFigure13b.png"
T2tt_1.figure  = "Aux.Fig.10b,Aux.Fig.13b"
#----exclusion source----
T2tt_1.addSource( 'obsExclusion', 'orig/T2tt_ATLAS_SUSY_2015_02_Obs.txt', 'txt', objectName = None, index = None )
T2tt_1.addSource( 'expExclusion', 'orig/T2tt_ATLAS_SUSY_2015_02_Exp.txt', 'txt', objectName = None, index = None )
T2tt_1.addSource( 'efficiencyMap', 'orig/EffMap_T2tt_SR2.txt', 'txt', objectName = None, index = None )
#T2bb_1.expExclusionM1.setSource( path, type, objectName = None, index = None )
#T2bb_1.expExclusionP1.setSource( path, type, objectName = None, index = None )
T2tt_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1469069/all"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR3")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR3", observedN = 1, expectedBG = 1.03 , bgError = 0.18)
T2tt = dataset.addTxName('T2tt')
T2tt.checked = 'NO'
T2tt.constraint ="[[['t']],[['t']]]"
T2tt.conditionDescription = None
T2tt.condition = None
T2tt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane( [[x,y]]*2 )
#---- new efficiency map -----
#----figure----
T2tt_1.figureUrl  = "http://hepdata.cedar.ac.uk/resource/1469069/figAuxiliaryFigure11b.png; http://hepdata.cedar.ac.uk/resource/1469069/figAuxiliaryFigure14b.png"
T2tt_1.figure  = "Aux.Fig.11b,Aux.Fig.14b"
#----exclusion source----
T2tt_1.addSource( 'obsExclusion', 'orig/T2tt_ATLAS_SUSY_2015_02_Obs.txt', 'txt', objectName = None, index = None )
T2tt_1.addSource( 'expExclusion', 'orig/T2tt_ATLAS_SUSY_2015_02_Exp.txt', 'txt', objectName = None, index = None )
T2tt_1.addSource( 'efficiencyMap', 'orig/EffMap_T2tt_SR3.txt', 'txt', objectName = None, index = None )
#T2bb_1.expExclusionM1.setSource( path, type, objectName = None, index = None )
#T2bb_1.expExclusionP1.setSource( path, type, objectName = None, index = None )
T2tt_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1469069/all"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("SR1")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "SR1", observedN = 12, expectedBG = 5.5 , bgError = 0.72)
T2tt = dataset.addTxName('T2tt')
T2tt.checked = 'NO'
T2tt.constraint ="[[['t']],[['t']]]"
T2tt.conditionDescription = None
T2tt.condition = None
T2tt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T2tt_1 = T2tt.addMassPlane( [[x,y]]*2 )
#---- new efficiency map -----
#----figure----
T2tt_1.figureUrl  = "http://hepdata.cedar.ac.uk/resource/1469069/figAuxiliaryFigure9b.png; http://hepdata.cedar.ac.uk/resource/1469069/figAuxiliaryFigure12b.png"
T2tt_1.figure  = "Aux.Fig.9b , Aux.Fig.12b"
#----exclusion source----
T2tt_1.addSource( 'obsExclusion', 'orig/T2tt_ATLAS_SUSY_2015_02_Obs.txt', 'txt', objectName = None, index = None )
T2tt_1.addSource( 'expExclusion', 'orig/T2tt_ATLAS_SUSY_2015_02_Exp.txt', 'txt', objectName = None, index = None )
T2tt_1.addSource( 'efficiencyMap', 'orig/EffMap_T2tt_SR1.txt', 'txt', objectName = None, index = None )
#T2bb_1.expExclusionM1.setSource( path, type, objectName = None, index = None )
#T2bb_1.expExclusionP1.setSource( path, type, objectName = None, index = None )
T2tt_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1469069/all"
databaseCreator.create()
