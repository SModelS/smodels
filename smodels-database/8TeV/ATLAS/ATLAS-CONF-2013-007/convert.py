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
info = MetaInfoInput('ATLAS-CONF-2013-007')
info.comment = 'A technical problem has been found in the pseudo-experiments used to make Table 4 in the conference note of March 1. (Partly) superseded by ATLAS-SUSY-2013-09, fastlim maps contain more topologies for ATLAS-CONF-2013-007. Topologies included here are, T1, T1bbbb, T1bbbt, T1bbqq, T1bbtt, T1btbt, T1btqq, T1bttt, T1qqtt, T1tttt, T2, T2bb, T2bt, T2tt, T5bbbb, T5bbbt, T5btbt, T5tbtb, T5tbtt, T5tttt, TGQ, TGQbbq, TGQbtq, TGQQtt.'
info.sqrts = '8.0'
info.private = False
info.lumi = '20.7'
info.url = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-007/'
info.prettyName = '2 SS leptons + 0-3 b-jets + Etmiss'
info.implementedBy = 'MT'


#+++++++ dataset block ++++++++++++++
dataset = DataSetInput('data')
dataset.setInfo(dataType = 'upperLimit', dataId = None)

#+++++++ next txName block ++++++++++++++
T6ttWW = dataset.addTxName('T6ttWW')
T6ttWW.checked ='VM'
T6ttWW.constraint ="[[['t+'],['W-']],[['t-'],['W+']]]"
T6ttWW.conditionDescription ="None"
T6ttWW.condition ="None"
T6ttWW.massConstraint = None
T6ttWW.source = 'ATLAS'
T6ttWWoff = dataset.addTxName('T6ttWWoff')
T6ttWWoff.constraint = "3.5*([[['t+'],['l-','nu']],[['t-'],['jet','jet']]] + [[['t-'],['l+','nu']],[['t+'],['jet','jet']]])"
T6ttWWoff.conditionDescription = "[[['t'],['mu','nu']],[['t'],['jet','jet']]]>[[['t'],['e','nu']],[['t'],['jet','jet']]]"
T6ttWWoff.condition = "Cgtr([[['t'],['mu','nu']],[['t'],['jet','jet']]],[[['t'],['e','nu']],[['t'],['jet','jet']]])"
T6ttWWoff.massConstraint = [['dm >= 169.0', 'dm <= 76.0'], ['dm >= 169.0', 'dm <= 76.0']]
T6ttWWoff.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T6ttWWLSP060 = T6ttWW.addMassPlane(2*[[x, y, 60.0]])
T6ttWWLSP060.figure = 'Fig. 15a'
T6ttWWLSP060.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-007/fig_15a.png'
T6ttWWLSP060.dataUrl = 'Not defined'
T6ttWWLSP060.setSources(dataLabels= ['obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T6ttWWLSP060.txt', 'orig/exclusionm1_T6ttWWLSP060.txt', 'orig/exclusionp1_T6ttWWLSP060.txt', 'orig/T6ttWWLSP060.txt'],
                 dataFormats= ['svg', 'svg', 'svg', 'txt'],units= [None, None, None, 'fb'])
T6ttWWoff.addMassPlane(T6ttWWLSP060)
#+++++++ next mass plane block ++++++++++++++
T6ttWWx200 = T6ttWW.addMassPlane(2*[[x, y, y/2.]])
T6ttWWx200.figure = 'Fig. 15b'
T6ttWWx200.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-007/fig_15b.png'
T6ttWWx200.dataUrl = 'Not defined'
T6ttWWx200.setSources(dataLabels= ['obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T6ttWWx200.txt', 'orig/exclusionm1_T6ttWWx200.txt', 'orig/exclusionp1_T6ttWWx200.txt', 'orig/T6ttWWx200.txt'],
                 dataFormats= ['svg', 'svg', 'svg', 'txt'],units= [None, None, None, 'fb'])
T6ttWWoff.addMassPlane(T6ttWWx200)

#+++++++ next txName block ++++++++++++++
T1tttt = dataset.addTxName('T1tttt')
T1tttt.checked ="VM"
T1tttt.constraint ="[[['t+','t-']],[['t-','t+']]]"
T1tttt.conditionDescription ="None"
T1tttt.condition ="None"
T1tttt.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1tttt = T1tttt.addMassPlane(2*[[x, y]])
T1tttt.figure = 'Fig. 8a'
T1tttt.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-007/fig_08a.png'
T1tttt.dataUrl = 'Not defined'
T1tttt.setSources(dataLabels= ['obsExclusion', 'obsExclusionM1', 'obsExclusionP1', 'upperLimits'],
                 dataFiles= ['orig/exclusion_T1tttt.txt', 'orig/exclusionm1_T1tttt.txt', 'orig/exclusionp1_T1tttt.txt', 'orig/T1tttt.txt'],
                 dataFormats= ['svg', 'svg', 'svg', 'txt'],units= [None, None, None, 'fb'])

#+++++++ next txName block ++++++++++++++
T5tttt = dataset.addTxName('T5tttt')
T5tttt.checked ="VM"
T5tttt.constraint ="[[['t+'],['t-']],[['t-'],['t+']]]+[[['t-'],['t+']],[['t-'],['t+']]]+[[['t+'],['t-']],[['t+'],['t-']]]"
T5tttt.conditionDescription ="None"
T5tttt.condition ="None"
T5tttt.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T5ttttLSP060 = T5tttt.addMassPlane(2*[[x, y, 60.0]])
T5ttttLSP060.figure = 'Fig. 9'
T5ttttLSP060.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-007/fig_09.png'
T5ttttLSP060.dataUrl = 'Not defined'
T5ttttLSP060.addSource(dataLabel= 'obsExclusion', dataFile = 'orig/exclusion_T5ttttLSP060.txt', dataFormat = 'svg')
T5ttttLSP060.addSource(dataLabel= 'obsExclusionM1', dataFile = 'orig/exclusionm1_T5ttttLSP060.txt', dataFormat = 'svg')
T5ttttLSP060.addSource(dataLabel= 'obsExclusionP1', dataFile = 'orig/exclusionp1_T5ttttLSP060.txt', dataFormat = 'svg')
T5ttttLSP060.addSource(dataLabel= 'upperLimits', dataFile = 'orig/T5ttttLSP060.txt', dataFormat = 'txt', unit = 'fb')

#+++++++ next txName block ++++++++++++++
T1btbt = dataset.addTxName('T1btbt')
T1btbt.checked ="VM"
T1btbt.constraint ="2*([[['t+','b']],[['t+','b']]]+[[['t-','b']],[['t-','b']]])"
T1btbt.conditionDescription ="None"
T1btbt.condition ="None"
T1btbt.source = 'ATLAS'
#+++++++ next mass plane block ++++++++++++++
T1btbt = T1btbt.addMassPlane(2*[[x, y]])
T1btbt.figure = 'Fig. 11'
T1btbt.figureUrl = 'https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2013-007/fig_11.png'
T1btbt.dataUrl = 'Not defined'
T1btbt.addSource(dataLabel= 'obsExclusion', dataFile = 'orig/exclusion_T1tbtb.txt', dataFormat = 'svg',
                    coordinateMap = {x : 0, y : 1, 'value': None})
T1btbt.addSource(dataLabel= 'obsExclusionM1', dataFile = 'orig/exclusionm1_T1tbtb.txt', dataFormat = 'svg')
T1btbt.addSource(dataLabel= 'obsExclusionP1', dataFile = 'orig/exclusionp1_T1tbtb.txt', dataFormat = 'svg')
T1btbt.addSource(dataLabel= 'upperLimits', dataFile = 'orig/T1tbtb.txt', dataFormat = 'txt', unit = 'pb', 
                  coordinateMap = {x : 0, y : 1, 'ul': 2}, scale = 0.001)




databaseCreator.create()
