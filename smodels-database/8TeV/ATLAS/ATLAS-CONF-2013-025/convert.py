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
info = MetaInfoInput('ATLAS-CONF-2013-025')
info.sqrts = '8.0'
info.private = False
info.lumi = '20.7'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-025/'
info.supersededBy = 'ATLAS-SUSY-2013-08' 
info.prettyName = '>= 5 (>=1 b-)  jets + 2,3 SF OS leptons + Etmiss'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T6ZZtt = dataset.addTxName('T6ZZtt')
T6ZZtt.checked ="AL"
T6ZZtt.constraint ="[[['Z'],['t']],[['Z'],['t']]]"
T6ZZtt.conditionDescription ="None"
T6ZZtt.condition ="None"
T6ZZtt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T6ZZttD180 = T6ZZtt.addMassPlane(2*[[x, y+180.0, y]])
T6ZZttD180.figure = 'Figure 10'
T6ZZttD180.figureUrl ='https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-025/figaux_10.png'
T6ZZttD180.dataUrl = 'Not defined'
T6ZZttD180.setSources(dataLabels= ['obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T6ZZttD180.txt', 'orig/exclusionm1_T6ZZttD180.txt', 'orig/exclusionp1_T6ZZttD180.txt', 'orig/T6ZZttD180.txt'],
                 dataFormats= ['svg', 'svg', 'svg', 'txt'],units= [None, None, None, 'fb'])



databaseCreator.create()
