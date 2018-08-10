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
info = MetaInfoInput('ATLAS-CONF-2013-048')
info.comment = 'T6bbWWM1300 can also be found in ATLAS-CONF-2013-065 as combination'
info.sqrts = '8.0'
info.private = False
info.lumi = '20.3'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-048/ http://cds.cern.ch/record/1547564'
info.supersededBy = 'ATLAS-SUSY-2013-19'
info.prettyName = '2 leptons + (b)jets + Etmiss'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T2bbWW = dataset.addTxName('T2bbWW')
T2bbWW.checked ="VM"
T2bbWW.constraint ="[[['b','W+']],[['b','W-']]]"
T2bbWW.conditionDescription ="None"
T2bbWW.condition ="None"
T2bbWW.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T2bbWW = T2bbWW.addMassPlane(2*[[x, y]])
T2bbWW.figure = 'Fig. 13f'
T2bbWW.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-048/fig_13f.png'
T2bbWW.dataUrl = 'Not defined'
T2bbWW.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T2bbWW.txt', 'orig/T2bbWW.dat'],
                 dataFormats= ['txt', 'txt'])

#+++++++ next txName block ++++++++++++++
T6bbWW = dataset.addTxName('T6bbWW')
T6bbWW.checked ="VM"
T6bbWW.constraint ="[[['b'],['W+']],[['b'],['W-']]]"
T6bbWW.conditionDescription ="None"
T6bbWW.condition ="None"
T6bbWW.source = "ATLAS"
T6bbWW.massConstraint = None
T6bbWWoff = dataset.addTxName('T6bbWWoff')
T6bbWWoff.constraint = "9.18*([[['b'],['L-','nu']],[['b'],['L+','nu']]])"
T6bbWWoff.conditionDescription = "None"
T6bbWWoff.condition = "None"
T6bbWWoff.massConstraint = [['dm >= 0.0', 'dm <= 76.0'], ['dm >= 0.0', 'dm <= 76.0']]
T6bbWWoff.source = "ATLAS"
#+++++++ next mass plane block ++++++++++++++
T6bbWWLSP001 = T6bbWW.addMassPlane(2*[[x, y, 1.0]])
T6bbWWLSP001.figure = 'Fig. 13a'
T6bbWWLSP001.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-048/fig_13a.png'
T6bbWWLSP001.dataUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-048/fig_13a_PRELIMINARY.data'
T6bbWWLSP001.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T6bbWWLSP001.txt', 'orig/T6bbWWLSP001.dat'],
                 dataFormats= ['txt', 'txt'])
T6bbWWoff.addMassPlane(T6bbWWLSP001)
#+++++++ next mass plane block ++++++++++++++
T6bbWWC150 = T6bbWW.addMassPlane(2*[[x, 150.0, y]])
T6bbWWC150.figure = 'Fig. 13e'
T6bbWWC150.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-048/fig_13e.png'
T6bbWWC150.dataUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-048/fig_13e_PRELIMINARY.data'
T6bbWWC150.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T6bbWWC150.txt', 'orig/T6bbWWC150.dat'],
                 dataFormats= ['txt', 'txt'])
T6bbWWoff.addMassPlane(T6bbWWC150)
#+++++++ next mass plane block ++++++++++++++
T6bbWWx200 = T6bbWW.addMassPlane(2*[[x, y*2.0, y]])
T6bbWWx200.figure = 'Fig. 13b'
T6bbWWx200.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-048/fig_13b.png'
T6bbWWx200.dataUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-048/fig_13b_PRELIMINARY.data'
T6bbWWx200.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T6bbWWx200.txt', 'orig/T6bbWWx200.dat'],
                 dataFormats= ['txt', 'txt'])
T6bbWWoff.addMassPlane(T6bbWWx200)
#+++++++ next mass plane block ++++++++++++++
T6bbWWD010 = T6bbWW.addMassPlane(2*[[x, x-10.0, y]])
T6bbWWD010.figure = "fig 13c"
T6bbWWD010.figureUrl = "https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-048/fig_13c.png"
T6bbWWD010.dataUrl = 'Not defined'
T6bbWWD010.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T6bbWWD010.txt', 'orig/T6bbWWD010.dat'],
                 dataFormats= ['txt', 'txt'])
T6bbWWoff.addMassPlane(T6bbWWD010)
#+++++++ next mass plane block ++++++++++++++
T6bbWWM1300 = T6bbWW.addMassPlane(2*[[300.0, x, y]])
T6bbWWM1300.figure = 'Fig. 13d'
T6bbWWM1300.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-048/fig_13d.png'
T6bbWWM1300.dataUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-048/fig_13d_PRELIMINARY.data'
T6bbWWM1300.setSources(dataLabels= ['obsExclusion', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T6bbWWM1300.txt', 'orig/T6bbWWM1300.dat'],
                 dataFormats= ['txt', 'txt'])
T6bbWWoff.addMassPlane(T6bbWWM1300)



databaseCreator.create()
