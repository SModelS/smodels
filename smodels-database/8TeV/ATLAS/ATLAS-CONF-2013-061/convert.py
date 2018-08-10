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
info = MetaInfoInput('ATLAS-CONF-2013-061')
info.comment = 'Superseeded by ATLAS-SUSY-2013-18 but publication has fewer topologies, e.g. T1btbt is onyl in conf note.'
info.sqrts = '8.0'
info.private = False
info.lumi = '20.1'
info.publication = 'http://cds.cern.ch/record/1557778'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-061/'
info.prettyName = '0 or >=1 leptons + jets + >= 3 b-jets + Etmiss'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T1bbbb = dataset.addTxName('T1bbbb')
T1bbbb.checked ="VM"
T1bbbb.constraint = "[[['b','b']],[['b','b']]]"
T1bbbb.conditionDescription ="None"
T1bbbb.condition ="None"
T1bbbb.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T1bbbb = T1bbbb.addMassPlane(2*[[x, y]])
T1bbbb.figure = 'Fig. 12a'
T1bbbb.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-061/fig_12a.png'
T1bbbb.dataUrl = 'Not defined'
T1bbbb.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T1bbbb.txt', 'orig/Fig12a.txt'],
                 dataFormats= ['txt', 'txt'],units= [None, 'fb'])

#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="VM"
T1tttt.constraint = "[[['t','t']],[['t','t']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane(2*[[x, y]])
T1tttt.figure = 'Fig. 12b'
T1tttt.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-061/fig_12b.png'
T1tttt.dataUrl = 'Not defined'
T1tttt.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T1tttt.txt', 'orig/Fig12b.txt'],
                 dataFormats= ['txt', 'txt'],units= [None, 'fb'])

#+++++++ next txName block ++++++++++++++
T1btbt = dataset.addTxName('T1btbt')
T1btbt.checked ="VM"
T1btbt.constraint = "[[['t','b']],[['t','b']]]"
T1btbt.conditionDescription ="None"
T1btbt.condition ="None"
T1btbt.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T1btbt = T1btbt.addMassPlane(2*[[x, y]])
T1btbt.figure = 'Fig. 12c'
T1btbt.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-061/fig_12c.png'
T1btbt.dataUrl = 'Not defined'
T1btbt.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T1tbtb.txt', 'orig/Fig12c.txt'],
                 dataFormats= ['txt', 'txt'],units= [None, 'fb'])



databaseCreator.create()
