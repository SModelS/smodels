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
info = MetaInfoInput('ATLAS-SUSY-2013-08')
info.comment = ""
info.sqrts = '8.0'
info.private = False
info.lumi = '20.3'
info.publication = 'http://link.springer.com/article/10.1140/epjc/s10052-014-2883-6'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-08/'
info.supersededBy = ""
info.arxiv = 'http://arxiv.org/abs/1403.5222'
info.contact = "ATLAS collaboration"
info.prettyName = 'Z + b-jets + Etmiss'
info.supersedes = 'ATLAS-CONF-2013-025'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T6ZZtt = dataset.addTxName('T6ZZtt')
T6ZZtt.constraint ="[[['Z'],['t']],[['Z'],['t']]]"
T6ZZtt.conditionDescription ="None"
T6ZZtt.condition ="None"
T6ZZtt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T6ZZttD180 = T6ZZtt.addMassPlane(2*[[x, y+180.0, y]])
T6ZZttD180.figure = "Fig 6"
T6ZZttD180.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/SUSY-2013-08/figaux_06.png"
T6ZZttD180.dataUrl = 'Not defined'
T6ZZttD180.setSources(dataLabels= ['expExclusion', 'expExclusionM1', 'expExclusionP1', 'obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/expectexclusion_T6ZZtt.txt', 'orig/expectexclusionm1_T6ZZtt.txt', 'orig/expectexclusionp1_T6ttZZ.txt', 'orig/exclusion_T6ZZtt.txt', 'orig/exclusionm1_T6ZZtt.txt', 'orig/exclusionp1_T6ZZtt.txt', 'orig/limit_T6ZZtt.txt'],
                 dataFormats= ['txt', 'txt', 'txt', 'txt', 'txt', 'txt', 'txt'],units= [None, None, None, None, None, None, 'fb'])



databaseCreator.create()
