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
info = MetaInfoInput('ATLAS-CONF-2012-105')
info.sqrts = '8.0'
info.private = False
info.lumi = '5.8'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2012-105/'
info.supersededBy = 'ATLAS-SUSY-2013-09'
info.prettyName = '2 SS leptons + >= 4 jets + Etmiss'
info.implementedBy = 'MT'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="AL"
T1tttt.constraint ="[[['t+','t-']],[['t+','t-']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane(2*[[x, y]])
T1tttt.figure = 'Fig. 3'
T1tttt.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2012-105/fig_03.png'
T1tttt.dataUrl = 'Not defined'
T1tttt.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/T1tttt_excl.txt', 'orig/T1tttt.txt'],
                 dataFormats= ['txt', 'txt'],units= [None, 'fb'])



databaseCreator.create()
