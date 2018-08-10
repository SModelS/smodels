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
info = MetaInfoInput('ATLAS-SUSY-2014-03')
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2014-03/'
info.sqrts = 8
info.lumi = 20.3 
info.prettyName = '>= 2(c-)jets + Etmiss'
info.arxiv = 'http://arxiv.org/abs/1501.01325'
info.publication = 'http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.114.161801'



#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("mCT200")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "mCT200", observedN = 11, expectedBG = 16 , bgError = 3, upperLimit = '3.760E-01*fb', expectedUpperLimit = '5.569E-01*fb')
#+++++++ next txName block ++++++++++++++
TScharm = dataset.addTxName('TScharm')
TScharm.checked = 'NO'
TScharm.constraint ="[[['c']],[['c']]]"
TScharm.conditionDescription = "None"
TScharm.condition ='None'
TScharm.source ='ATLAS'
#+++++++ next mass plane block ++++++++++++++
TScharm_1 = TScharm.addMassPlane([[x,y]]*2)
TScharm_1.figure  = "Figure9b"
TScharm_1.figureUrl  = "http://hepdata.cedar.ac.uk/resource/1337472/figAuxiliaryFigure9b.png"
TScharm_1.addSource('obsExclusion', 'orig/TScharm_Exclusion.txt', 'txt')
TScharm_1.addSource('efficiencyMap', 'orig/mCT200.dat', 'txt')
TScharm_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1337472/first"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("mCT150")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "mCT150", observedN = 19, expectedBG = 30 , bgError = 6, upperLimit = '5.018E-01*fb', expectedUpperLimit = '8.395E-01*fb')
#+++++++ next txName block ++++++++++++++
TScharm = dataset.addTxName('TScharm')
TScharm.checked = 'NO'
TScharm.constraint ="[[['c']],[['c']]]"
TScharm.conditionDescription = "None"
TScharm.condition ='None'
TScharm.source ='ATLAS'
#+++++++ next mass plane block ++++++++++++++
TScharm_1 = TScharm.addMassPlane([[x,y]]*2)
TScharm_1.figure  = "Figure2"
TScharm_1.figureUrl  = "http://hepdata.cedar.ac.uk/resource/1337472/figFigure2.png"
TScharm_1.addSource('obsExclusion', 'orig/TScharm_Exclusion.txt', 'txt')
TScharm_1.addSource('efficiencyMap', 'orig/mCT150.dat', 'txt')
TScharm_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1337472/first"


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput("mCT250")
dataset.setInfo(dataType = 'efficiencyMap', dataId = "mCT250", observedN = 4, expectedBG = 8.2 , bgError = 1.9, upperLimit = '2.544E-01*fb', expectedUpperLimit = '4.027E-01*fb')
#+++++++ next txName block ++++++++++++++
TScharm = dataset.addTxName('TScharm')
TScharm.checked = 'NO'
TScharm.constraint ="[[['c']],[['c']]]"
TScharm.conditionDescription = "None"
TScharm.condition ='None'
TScharm.source ='ATLAS'
#+++++++ next mass plane block ++++++++++++++
TScharm_1 = TScharm.addMassPlane([[x,y]]*2)
TScharm_1.figure  = "Figure9c"
TScharm_1.figureUrl  = "http://hepdata.cedar.ac.uk/resource/1337472/figAuxiliaryFigure9c.png"
TScharm_1.addSource('obsExclusion', 'orig/TScharm_Exclusion.txt', 'txt')
TScharm_1.addSource('efficiencyMap', 'orig/mCT250.dat', 'txt')
TScharm_1.dataUrl  = "http://hepdata.cedar.ac.uk/view/ins1337472/first"

databaseCreator.create()
